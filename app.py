import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from chembl_webresource_client.new_client import new_client

# ---------------------------
# Page setup
# ---------------------------
st.set_page_config(page_title="ChEMBL Substructure Search", layout="wide")
st.title("🔬 ChEMBL Substructure Search App")

# ---------------------------
# Cached target name fetcher
# ---------------------------
@st.cache_data
def get_target_name(tid):
    try:
        res = new_client.target.get(tid)
        return res.get("pref_name") or tid
    except:
        return tid

# ---------------------------
# SMARTS input
# ---------------------------
smarts = st.text_input("Enter SMARTS / substructure", "c1ccccc1")
query_mol = Chem.MolFromSmarts(smarts)

if smarts and query_mol is None:
    st.error("Invalid SMARTS string")
    st.stop()

# ---------------------------
# Target selection
# ---------------------------
st.subheader("🎯 Target Selection")

mode = st.radio(
    "Select mode",
    ["Manual (multi-ID)", "Search & auto-expand"]
)

targets = {}

# ---------------------------
# MODE 1: Manual multi-ID
# ---------------------------
if mode == "Manual (multi-ID)":

    target_text = st.text_area(
        "Enter targets and ChEMBL IDs (one per line)",
        """CCR5: CHEMBL3473, CHEMBL3217397
CCR3: CHEMBL240"""
    )

    for line in target_text.split("\n"):
        if ":" not in line:
            continue

        name, ids = line.split(":")
        name = name.strip()

        id_list = [i.strip() for i in ids.split(",") if i.strip()]

        if id_list:
            targets[name] = id_list

# ---------------------------
# MODE 2: Auto expand
# ---------------------------
else:
    target_client = new_client.target

    search_term = st.text_input("Search target (e.g. CCR5)")
    include_family = st.checkbox("Include PROTEIN FAMILY", value=True)

    if search_term:
        results = target_client.search(search_term)

        grouped = {}

        for r in results:
            if not r.get("target_chembl_id"):
                continue

            if r.get("organism") != "Homo sapiens":
                continue

            ttype = r.get("target_type")

            if ttype == "SINGLE PROTEIN" or (include_family and ttype == "PROTEIN FAMILY"):
                name = r.get("pref_name") or search_term

                if name not in grouped:
                    grouped[name] = []

                grouped[name].append(r["target_chembl_id"])

        options = {
            f"{k} ({len(v)} IDs)": (k, v)
            for k, v in grouped.items()
        }

        selected = st.multiselect("Select targets", list(options.keys()))

        for s in selected:
            name, ids = options[s]
            targets[name] = ids

# ---------------------------
# Filters
# ---------------------------
st.subheader("⚙️ Filters")

col1, col2 = st.columns(2)

with col1:
    activity_types = st.multiselect(
        "Activity types",
        ["IC50", "EC50"],
        default=["IC50", "EC50"]
    )

with col2:
    max_value = st.number_input(
        "Max activity (nM)",
        value=10000
    )

# ---------------------------
# SAR options
# ---------------------------
st.subheader("🧬 SAR Options")

use_pactivity = st.checkbox("Convert to pActivity (-log10 M)", value=True)

aggregation_method = st.selectbox(
    "Aggregation method",
    ["min", "mean", "median"]
)

# ---------------------------
# Run search
# ---------------------------
if st.button("🚀 Run Search"):

    if not targets:
        st.warning("Please select at least one target")
        st.stop()

    activity = new_client.activity
    results = []

    st.info("Fetching data from ChEMBL...")
    progress = st.progress(0)

    total_targets = len(targets)

    # ---------------------------
    # Loop through targets
    # ---------------------------
    for idx, (target_name, target_ids) in enumerate(targets.items()):

        st.write(f"🔎 Processing {target_name}")

        fetched_count = 0
        matched_count = 0

        for target_id in target_ids:

            st.write(f"   ↳ querying {target_id}")

            try:
                acts = activity.filter(
                    target_chembl_id=target_id,
                    standard_type__in=activity_types,
                    assay_type="B",
                    standard_value__lte=max_value
                ).only([
                    "molecule_chembl_id",
                    "canonical_smiles",
                    "standard_type",
                    "standard_value",
                    "standard_units"
                ])

                for i, a in enumerate(acts):
                    fetched_count += 1

                    if i % 500 == 0 and i > 0:
                        st.write(f"{target_name}: processed {i}")

                    smiles = a.get("canonical_smiles")
                    if not smiles:
                        continue

                    mol = Chem.MolFromSmiles(smiles)

                    if mol and mol.HasSubstructMatch(query_mol):

                        matched_count += 1

                        results.append({
                            "target_name": target_name,
                            "target_id": target_id,
                            "molecule_chembl_id": a.get("molecule_chembl_id"),
                            "smiles": smiles,
                            "standard_type": a.get("standard_type"),
                            "standard_value": a.get("standard_value"),
                            "standard_units": a.get("standard_units"),
                        })

            except Exception as e:
                st.warning(f"Error with {target_id}")
                continue

        st.write(f"✅ {target_name}: fetched {fetched_count} | matched {matched_count}")

        progress.progress((idx + 1) / total_targets)

    # ---------------------------
    # DataFrame
    # ---------------------------
    df = pd.DataFrame(results)

    if df.empty:
        st.warning("No matches found")
        st.stop()

    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")

    # Remove duplicates across multiple IDs
    df = df.drop_duplicates(
        subset=["molecule_chembl_id", "target_name", "standard_value"]
    )

    # ---------------------------
    # pActivity conversion
    # ---------------------------
    if use_pactivity:
        df["activity"] = -np.log10(df["standard_value"] * 1e-9)
    else:
        df["activity"] = df["standard_value"]

    # ---------------------------
    # Output: long format
    # ---------------------------
    st.success(f"✅ Total matches: {len(df)}")

    st.write("### Results per target")
    st.write(df.groupby("target_name").size())

    st.write("### 📋 Raw Data")
    st.dataframe(df, use_container_width=True)

    csv = df.to_csv(index=False).encode("utf-8")

    st.download_button(
        "⬇️ Download Raw CSV",
        csv,
        "chembl_results.csv",
        "text/csv"
    )

    # ---------------------------
    # SAR MATRIX
    # ---------------------------
    st.write("## 🧬 SAR Matrix")

    sar_df = df.pivot_table(
        index=["smiles", "molecule_chembl_id"],
        columns="target_name",
        values="activity",
        aggfunc=aggregation_method
    ).reset_index()

    sar_df.columns.name = None

    st.dataframe(sar_df, use_container_width=True)

    csv_sar = sar_df.to_csv(index=False).encode("utf-8")

    st.download_button(
        "⬇️ Download SAR Matrix",
        csv_sar,
        "sar_matrix.csv",
        "text/csv"
    )