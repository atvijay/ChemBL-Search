# ChemBL-Search
A Streamlit-based web application for performing substructure searches on ChEMBL bioactivity data, with filtering by target, activity type (IC50/EC50), and potency thresholds.

This tool enables medicinal chemists to quickly identify compounds containing a given substructure and analyze their activity against selected biological targets.

🚀 Features
Substructure search (SMARTS-based) using RDKit
🎯 Multi-target querying
Enter ChEMBL target IDs (e.g., PIM1, PIM3)
Or search targets by name
⚙️ Activity filtering
IC50 / EC50 selection
Custom activity cutoff (nM)
📊 Live progress tracking
Records fetched and matched per target
🧠 Automatic target name resolution
Converts ChEMBL IDs → protein names (cached for speed)
📁 Downloadable results (CSV)
📈 Per-target result summary

The app will return all matching compounds with associated activity data.

🛠️ Tech Stack
Python
Streamlit (UI)
RDKit (substructure matching)
ChEMBL Web Resource Client (data access)
Pandas (data handling)
📦 Installation

Clone the repository:

git clone https://github.com/your-username/chembl-substructure-search.git
cd chembl-substructure-search

Create environment (recommended):

conda create -n chembl_env python=3.10
conda activate chembl_env

Install dependencies:

pip install -r requirements.txt
▶️ Running the App
streamlit run app.py

The app will open in your browser (usually http://localhost:8501).

🧾 Input Guide
1. Substructure
Enter a SMILES string


2. Target Selection

Option A: Enter ChEMBL IDs

CHEMBL1111, CHEMBL2222
Option B: Search by name and select from list
3. Filters
Activity type: IC50 / EC50
Max activity value (nM)
📊 Output

The app returns a table with:

Column	Description
target_name	Protein name 
target_id	ChEMBL target ID
molecule_chembl_id	Compound ID
smiles	Molecular structure
standard_type	IC50 / EC50
standard_value	Activity value
standard_units	Typically nM

Results can be downloaded as a CSV file.

⚡ Performance Notes
Large targets may take several minutes (thousands of compounds)
Progress is shown during execution
Target names are cached to reduce API calls
⚠️ Known Limitations
Depends on ChEMBL API availability (may occasionally return 500 errors)
Substructure matching is done locally (can be slow for large datasets)
Only binding assays (assay_type = "B") are included
🔮 Future Improvements
Add pIC50 / pEC50 calculation
Include UniProt ID / gene symbols
Structure visualization (RDKit images)
Similarity search (Tanimoto)
Scaffold clustering / SAR analysis
Local caching of ChEMBL data
🤝 Contributing

Contributions are welcome!
Feel free to open issues or submit pull requests.

📜 License

MIT License

🙌 Acknowledgements
ChEMBL
RDKit
Streamlit

👤 Author

Vijayendar Yedulla

Medicinal Chemist 
