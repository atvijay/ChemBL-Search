import streamlit as st
import pandas as pd
from rdkit import Chem
from chembl_webresource_client.new_client import new_client

# ---------------------------
# Page setup
# ---------------------------
st.set_page_config(page_title="ChEMBL Substructure Search", layout="wide")

st.title("🔬 ChEMBL Substructure Search App")

# ---------------------------
# Input: SMARTS
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

mode = st.radio("Select mode", ["By ChEMBL ID", "Search by name"])

targets = {}

if mode == "By ChEMBL ID":
    target_input = st.text_input(
        "Enter target ChEMBL IDs (comma separated)",
        "CHEMBL3473, CHEMBL3217397"
    )

    # ✅ FIXED: unique readable labels
    targets = {
        f"Target_{i+1} ({tid.strip()})": tid.strip()
        for i, tid in enumerate(target_input.split(","))
        if tid.strip()
    }

else:
    target_client = new_client.target

    search_term = st.text_input("Search target (e.g. CCR5)")

    if search_term:
        results = target_client.search(search_term)

        options = {
            f"{r['pref_name']} ({r['target_chembl_id']})": r['target_chembl_id']
            for r in results
            if r.get("target_chembl_id")
            and r.get("target_type") == "SINGLE PROTEIN"
            and r.get("organism") == "Homo sapiens"
        }

        selected = st.multiselect("Select targets", list(options.keys()))
        targets = {s: options[s] for s in selected}

# ---------------------------
# Query settings
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
# Run button
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
    for idx, (target_name, target_id) in enumerate(targets.items()):

        st.write(f"🔎 Processing {target_name}")

        fetched_count = 0
        matched_count = 0

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
                    st.write(f"{target_name}: processed {i} records")

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
            st.error(f"❌ Error with {target_name}: {str(e)}")
            continue

        # ✅ Debug summary per target
        st.write(
            f"✅ {target_name}: fetched {fetched_count} records | matched {matched_count}"
        )

        progress.progress((idx + 1) / total_targets)

    # ---------------------------
    # DataFrame
    # ---------------------------
    df = pd.DataFrame(results)

    # Optional sorting
    if not df.empty:
        df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
        df = df.sort_values("standard_value")

    # ---------------------------
    # Output
    # ---------------------------
    st.success(f"✅ Total matches: {len(df)}")

    # ✅ Show distribution by target
    if not df.empty:
        st.write("### Results per target")
        st.write(df.groupby("target_name").size())

        st.dataframe(df, use_container_width=True)

        csv = df.to_csv(index=False).encode("utf-8")

        st.download_button(
            "⬇️ Download CSV",
            csv,
            "chembl_results.csv",
            "text/csv"
        )
    else:
        st.warning("No matches found")