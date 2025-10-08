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
    parser.add_argument('--target', type=str, required=True, help='The targeting compound supported: 2,3-BDO, fatty_acid')
    parser.add_argument('--solvent', type=str, required=True, help='The solvent to be predicted')
    args = parser.parse_args()
    
    # Check input
    if Chem.MolFromSmiles(args.solvent): # if input is a valid SMILES
        solvent_smiles = args.solvent
    else:
        try:
            solvent_smiles = pcp.get_compounds(args.solvent, 'name')
            solvent_smiles = solvent_smiles[0].isomeric_smiles
        except ValueError as e:
            print("Input cannot be covnerted to SMILES. Error message: ", e)
    
    # Compute result
    if args.target =='2,3-BDO':
        model_path = 'Diol/model_bdo_manually.pkl'
        scale_path = 'Diol/scaler_bdo_manually.pkl'
        featu_path = 'Diol/feature_bdo_manually.json'
        
        model_bdo = joblib.load(model_path)
        scaler_bdo = joblib.load(scale_path)
        with open(featu_path, 'r') as f1:
            feature_bdo = json.load(f1)
        
        temp_K = 273.15+25 
        fps_bdo, sim_bdo = predict_rdkfp( list(solvent_smiles), feature_bdo, model_bdo, scaler=scaler_bdo)
        sim = -sim_bdo/(2.303*1.987*0.001*temp_K)  # logP
        sim = np.power( 10, sim)  #sim = np.exp(sim)
        
        fitted_a, fitted_b = 0.057, 0.021
        sim = fitted_a * sim + fitted_b
        
        print( 'Predicted result: ', sim )
    if args.target =='fatty_acid':
        print('Under construction')
    else:
        print( f"Target compound {args.target} is not supported yet...")
        
        
if __name__ == "__main__":
    main()