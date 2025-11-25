# RRAM_probabilistic_computer_for_molecular_docking

**MATLAB codes for coordinate extraction, distance computation, QUBO/adjacency construction, Gaussian CDF/sigmoid fitting, and probabilistic MWCP verification** used in the study:

**“A Hardware Demonstration of a Universal Programmable RRAM-based Probabilistic Computer for Molecular Docking”.**

This repository contains illustrative MATLAB scripts supporting the coordinate extraction, distance computation, QUBO formulation, and verification of the 42-node Maximum Weighted Clique Problem (MWCP) derived from a molecular docking instance (LolA–LolCDE complex, PDB ID: *7arm*).
These scripts provide a transparent demonstration of the computational workflow.  


---

## 📁 Repository Structure
```

├── README.md # Repository documentation
├── 7arm.pdb # PDB file used in this study
├── Extract_coordinates_pairwise_distance.m   # Extract pharmacophore 3D coordinates and distance matrix
├── Reorder_coordinates_pairwise_distance.m   # Standardize pharmacophore ordering and recompute distances 
├── Compatibility_check_adjacency_matrix.m    # Build adjacency/QUBO matrix from pharmacophore point binding pairs
├── MWCP_solver_42nodes_pcomputing.m          # Probabilistic 42-node MWCP solver (verification for maximum weighted clique solutions)
└── GaussianCDF_vs_Sigmoid_Fitting.m          # Gaussian CDF vs sigmoid fitting 

```
## 🚀 Running the Code
Install MATLAB R2020+ or later

No additional toolboxes required

Scripts tested on Windows 11

Run the following scripts in order (you may skip steps depending on your purpose):

### **1. Coordinate extraction and pairwise distances**

```matlab
Extract_coordinates_pairwise_distance.m
Reorder_coordinates_pairwise_distance.m
```
### 2. Adjacency/QUBO matrix construction
```matlab
Compatibility_check_adjacency_matrix.m
```
### 3. MWCP probabilistic solver (verifying MWC solutions of the 42-node MWCP)
```matlab
MWCP_solver_42nodes_pcomputing.m
```
### 4. Visualization of fittings between Gaussian CDFs and sigmoid curves 
```matlab
GaussianCDF_vs_Sigmoid_Fitting.m
```


## 📜 License
This project is released under the MIT License.
You are free to use, modify, and distribute the code with proper attribution.

## 📬 Contact
For questions regarding the code, please contact:
Yihan He, email: [yihan_he@u.nus.edu]

## 📚 How to Cite
If you use this code in your research, please cite:

Code repository:
RRAM_probabilistic_computer_for_molecular_docking, GitHub, 2025.

DOI: (to be added after Zenodo release)
