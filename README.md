# LEAPS

LEAPS: Liquid–Liquid Extraction Artificial Predictor and Screener for biomass separation

This project aims to develop fast-screening tools for identifying solvent candidates in the liquid-liquid extraction of biomass downstream separation. 

An overall workflow:
![An overall of the workflow](./TOC.png)

Example usage: 
```bash
python Predictor.py --target 2,3-butanediol --solvent 1-butanol 1-hexanol 2-methyl-1-butanol 2-ethyl-1-pentanol
```
returns:
```bash
Final input SMILES representations are:  ['CCCCO', 'CCCCCCO', 'CCC(C)CO', 'CCCC(CC)CO']
========= Predicted result (Kd, Input) ========
0.756  1-butanol
0.406  1-hexanol
0.686  2-methyl-1-butanol
0.272  2-ethyl-1-pentanol
```

Current compounds with experimental validations:

1. Diol: 2,3-butanediol (2,3-BDO). Related publications: [1] Ind. Eng. Chem. Res. 2025, 64, 32, 15790–15799 (https://doi.org/10.1021/acs.iecr.5c01569) [2] Chem Bio Eng. 2025, 2, 4, 210–228 (https://pubs.acs.org/doi/10.1021/cbe.4c00170)

   Details:
       - Code_ML_screening.ipynb: Codes containing data processing, ML parameterization, solvent screening, and visualization
       
       - CSV files: The DFT comptuted data and experimental data for validation.
       
       - JSON and PKL files: The parameterized model and parameters ready for prediction

2. Long-chain fatty acid: under investigation.

This project is supported by U.S. Department of Energy, Bioenergy Technologies Office (BETO), via the Bioprocessing Separation Consortium. 
