
# BindCraft Output Analysis Web App (BOA)

This is a Streamlit-based web application for analyzing **BindCraft** output, exploring **interface metrics**, and visualizing **3D protein structures**.

This tool helps you interpret the results generated after running BindCraft, including the Accepted PDB files and the `final_design_stats.csv` file.

For more information about the output parameters of BindCraft, this link is so useful: 
https://yarrowmadrona.medium.com/denovo-protein-binder-design-using-bindcraft-part-ii-you-have-your-data-now-what-db9288e782b3


---

## Getting Started
### Accessing the bindcraft-output-analysis web app
This application is hosted on Streamlit Cloud, providing easy access without any local installation required. The link :
https://bindcraftoutputanalyse-ednrdmyx97eon4hlfkjatb.streamlit.app/

### Generating BindCraft Analysis Output
To utilize this analysis app, you first need to generate output from BindCraft. The original BindCraft software repository can be found here: https://github.com/martinpacesa/BindCraft.


## 🚀 Features

### ✅ 1. Summary Tab
- Displays all analyzed results in a single interactive table  
- Combines your original CSV with additional metrics computed by MDAnalysis  
- Provides buttons to download:
  - Updated summary CSV  
  - Ranked designs CSV  

---

### ✅ 2. Visualization Tab
- Interactive **3D structure viewer** using *py3Dmol*
- Customizable chain coloring  
- Counts and plots the number of **4 Å contacts** between binder and target  

---

### ✅ 3. Details Tab
For each design, you can inspect:
- Target and binder sequences  
- Residue–residue contact pairs (3 Å and 4 Å)  
- H-bond–like interactions  
- Interface residue counts   

---

## 📁 Input Requirements

Upload the following:

1. **Accepted PDB designs**  
   Output from BindCraft.

2. **`final_design_stats.csv`**  
   The BindCraft metrics file.

The app will automatically match each PDB to its corresponding row in the CSV.

---

## 🛠️ How to use it?

You have two choices:  
- install all the dependencies.
- just use the Wep App via the link.  (https://bindcraftoutputanalyse-ednrdmyx97eon4hlfkjatb.streamlit.app/)


If you want to install it locally:
####  Dependencies 

git clone https://github.com/mj72git/BindCraft_Output_Analyse.git
cd BindCraft_Output_Analyse-main

- Python
- MDAnalysis
- numpy
- pandas
- plotly
- py3Dmol
- streamlit




## References
[1] Pacesa, M. (2024). BindCraft: one-shot design of functional protein binders. bioRxiv. [https://www.biorxiv.org/content/10.1101/2024.09.30.615802v1]

[2] A-Yarow github page: https://github.com/A-Yarrow
