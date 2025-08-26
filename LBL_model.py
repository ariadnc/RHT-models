"""
Created on Tue Jun 24 01:43:56 2025

Compute emissivity profile of a single file with the LBL model

@author: ariad
"""

import pandas as pd
import numpy as np
from radis import calc_spectrum
from scipy.constants import pi, sigma as Sigma  #Stefan-Boltzmann constant

def compute_radiative_heat_transfer_rte(file_path,
                                        sheet_name=0,
                                        distance_col="Distance (cm)",
                                        temperature_col="Temperature (K)",
                                        fv_col="Particle_volume_fraction (cm3/cm3)",
                                        molecules=['OH', "H2O", 'CO', 'CO2'],
                                        spectral_range=(300, 10000),
                                        database='hitemp',
                                        C_eta=1182,
                                        C_par=7.0,
                                        path_length=1.0):   # cm
    
    df = pd.read_excel(file_path, sheet_name=sheet_name)
    T_profile = df[temperature_col].values
    f_v_profile = df[fv_col].values if fv_col in df.columns else np.zeros(len(df))

    mole_fractions = {}
    
    for mol in molecules:
        col = f"Mole_fraction_{mol} ()"
        if col not in df.columns:
            raise ValueError(f"Missing column for {mol}")
        mole_fractions[mol] = df[col].values

    x_vals = df[distance_col].values
    emissivity_list = []

    for i in range(len(x_vals)):
        print(f'We are at this point {i + 1}/{len(x_vals)}')
        
        T = T_profile[i]
        f_v = f_v_profile[i]
        I_L_total = 0

        for mol in molecules:
            mf = mole_fractions[mol][i]
            if mf <= 0:
                continue

            spec = calc_spectrum(spectral_range[0], spectral_range[1],
                                 molecule=mol,
                                 databank=database,
                                 Tgas=T,
                                 isotope=1,
                                 path_length=path_length,
                                 mole_fraction=mf,
                                 pressure=1, #atm
                                 wstep='auto',
                                 cutoff=1e-30)
          
            I_L_total += spec.get_integral("radiance_noslit")

        emissivity_gas = (I_L_total * pi) / (Sigma * T**4)

        k_soot = C_par * f_v * T
        emissivity_soot = 1 - np.exp(-k_soot * path_length*0.01)
      
        emissivity_total = emissivity_gas + emissivity_soot
        emissivity_list.append(emissivity_total)
        
    df["Emissivity (-)"] = emissivity_list

    return df, x_vals, np.array(emissivity_list)

if __name__ == "__main__":
    df, x_vals, emissivity_array = compute_radiative_heat_transfer_rte(file_path=r'C:\Users\ariad\OneDrive\Documentos\MATLAB\H2_1atm\excelCO2\H2_1atm_40_O2_60_CO2.xlsx')