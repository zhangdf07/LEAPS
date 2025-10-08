# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 10:23:29 2025

@author: d2j
"""

import argparse, json, joblib
import pandas as pd
import numpy as np

import pubchempy as pcp

from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors


import warnings
warnings.filterwarnings("ignore")


def predict_rdkfp(smi, fp_name, model, scaler=None, fp_calc=None):
    """
    smi = a list of SMILES
    """
    if fp_calc==None:
        fp_calc = MoleculeDescriptors.MolecularDescriptorCalculator(fp_name)
    fps = [list( fp_calc.CalcDescriptors( Chem.MolFromSmiles(s) ) ) for s in smi]
    new_df = pd.DataFrame( data=fps, columns=fp_name )
    if scaler is not None:
        new_df = scaler.transform(new_df)
    try:
        predict_new = model.predict(new_df)
    except:
        print('Prediction issue for', smi, new_df)
        predict_new = None
    
    return new_df, predict_new


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--target', type=str, default="2,3-butanediol", help='The targeting compound supported: 2,3-BDO, fatty_acid')
    parser.add_argument('--solvent', nargs='*', help='The solvents to be predicted')
    args = parser.parse_args()
    
    #print('Original input is: ', args.solvent )
    
    # Check input
    input_raw, input_smi, input_failed = [],[],[]
    for sol in args.solvent:
        try:
            solvent_smiles = pcp.get_compounds(sol, 'name')
            solvent_smiles = solvent_smiles[0].smiles # Older version: isomeric_smiles
            #print("Solvent input is used as a valid solvent name")
            input_raw.append(sol)
            input_smi.append(solvent_smiles)
        except:
            if Chem.MolFromSmiles(sol):
                solvent_smiles = args.solvent
                #print("Solvent input is used as a valid SMILES")
                input_raw.append(sol)
                input_smi.append(solvent_smiles)
            else:
                #raise ValueError("Input cannot be covnerted to SMILES.")
                input_failed.append(sol)
    print('Final valid input SMILES: ', input_smi )
    
    # Compute result
    if args.target =='2,3-butanediol':
        model_path = 'Diol/model_bdo_manually.pkl'
        scale_path = 'Diol/scaler_bdo_manually.pkl'
        featu_path = 'Diol/feature_bdo_manually.json'
        
        model_bdo = joblib.load(model_path)
        scaler_bdo = joblib.load(scale_path)
        with open(featu_path, 'r') as f1:
            feature_bdo = json.load(f1)
        
        temp_K = 273.15+25 
        fps_bdo, sim_bdo = predict_rdkfp( input_smi, feature_bdo, model_bdo, scaler=scaler_bdo)
        sim = -sim_bdo/(2.303*1.987*0.001*temp_K)  # logP
        sim = np.power( 10, sim)  #sim = np.exp(sim)
        
        fitted_a, fitted_b = 0.057, 0.021
        sim = np.round( fitted_a * sim + fitted_b, 3 )
        
        print( '========= Predicted results (Kd, Input) ========' )
        for raw, smi, result, in zip(input_raw, input_smi, sim):
            print( str(result).ljust(6), raw )
        #print( '========= Predicted result (Kd, Input) ========' )
        
    elif args.target =='fatty_acid':
        print('Under construction')
        
    else:
        print( f"Target compound {args.target} is not supported yet...")
        
        
if __name__ == "__main__":
    main()