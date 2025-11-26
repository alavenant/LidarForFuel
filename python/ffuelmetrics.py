"""
Fuel metrics LiDAR

This module provides functions to compute fuel metrics from bulk density profiles.
Produces similar outputs to the fCBDprofile_fuelmetrics function. This function is
intended for use with bulk density profiles derived from field data, MLS, TLS, or model outputs.
"""

import numpy as np
import pandas as pd
import warnings
from typing import Union


def ffuelmetrics(profile_table, threshold=0.012, H_PAI=0):
    """
    Compute fuel metrics from a bulk density profile.
    
    Parameters
    ----------
    profile_table : pandas.DataFrame
        A table with at least two columns: H (height of the strata) and BD (bulk density in kg/m³).
        If BD is not available, provide three more columns: PAD (plant area density in m²/m³),
        LMA (leaf mass area in g/m², use average value 141 if unknown), and WD (wood density
        in kg/m³, use average value 591 if unknown)
    threshold : numeric or str, default 0.012
        A critical bulk density threshold used to identify different strata limits such as
        midstorey height, canopy base, and canopy top. Can be either:
        - Numeric: a bulk density value (in kg/m³)
        - Str: a percentage of the maximum CBD value (e.g., "10%")
    H_PAI : numeric, default 0
        Height from which PAI and VCI_PAD is estimated. Default is 0 (means from ground to top).
        
    Returns
    -------
    tuple
        A tuple of two elements:
        1) A dict with all fuel metrics
        2) A pandas DataFrame with the PAD and CBD profile values plus additional columns
           if BD not provided (i.e., BD_rollmean, LMA, WD, WMA, FMA if BD not provided)
    """
    PAD_CBD_Profile = profile_table.copy()
    
    # Check if BD column exists
    if 'BD' not in PAD_CBD_Profile.columns and 'CBD' not in PAD_CBD_Profile.columns:
        # Need to compute BD from PAD, LMA, and WD
        if 'PAD' not in PAD_CBD_Profile.columns or 'LMA' not in PAD_CBD_Profile.columns or 'WD' not in PAD_CBD_Profile.columns:
            raise ValueError(
                "No BD (bulk density) column found. Please provide LMA and WD columns "
                "in order to compute BD. Use value LMA = 141 and WD = 591 if unknown"
            )
        
        # LMA g/m² => kg/m²
        PAD_CBD_Profile['LMA'] = PAD_CBD_Profile['LMA'] / 1000
        
        # Surface volume ratio (SVR: m²/m³) for 4mm diameter twigs
        SVR = 1 / 0.002
        # Wood mass area (WMA)
        PAD_CBD_Profile['WMA'] = PAD_CBD_Profile['WD'] / SVR
        
        # Partition of wood and leaves (M. Soma phd thesis data)
        partW = 0.51
        partL = 0.49
        
        # Fuel mass area
        PAD_CBD_Profile['FMA'] = 1 / ((partW / PAD_CBD_Profile['WMA']) + 
                                       (partL / PAD_CBD_Profile['LMA']))
        
        # CBD in kg/m³
        PAD_CBD_Profile['CBD'] = PAD_CBD_Profile['PAD'] * PAD_CBD_Profile['FMA']
    elif 'BD' in PAD_CBD_Profile.columns:
        # Rename BD to CBD for consistency
        PAD_CBD_Profile['CBD'] = PAD_CBD_Profile['BD']
    
    # Get strata depth
    if len(PAD_CBD_Profile) > 1:
        d = PAD_CBD_Profile['H'].iloc[1] - PAD_CBD_Profile['H'].iloc[0]
    else:
        d = 0.5  # Default value
    
    # Define threshold when threshold is a proportion of CBD max
    if isinstance(threshold, str) and '%' in threshold:
        threshold_prop = float(threshold.replace('%', '')) / 100
        threshold = np.max(PAD_CBD_Profile[PAD_CBD_Profile['H'] > 1]['CBD']) * threshold_prop
    else:
        threshold = float(threshold)
    
    # No data above 0.5m
    if np.max(PAD_CBD_Profile['H']) < 0.5:
        warnings.warn("NULL return: no data above 0.5m height")
        return _create_null_return_ffuelmetrics(PAD_CBD_Profile)
    
    # Get CBD roll mean to smooth the profiles & get the profile above CBD threshold
    if len(PAD_CBD_Profile) > 3:
        PAD_CBD_Profile['CBD_rollM'] = PAD_CBD_Profile['CBD'].rolling(window=3, center=False).mean()
        PAD_CBD_Profile['CBD_rollM'].iloc[:3] = PAD_CBD_Profile['CBD'].iloc[:3]
        PAD_CBD_Profile_threshold = PAD_CBD_Profile[PAD_CBD_Profile['CBD_rollM'] > threshold].copy()
    else:
        PAD_CBD_Profile['CBD_rollM'] = PAD_CBD_Profile['CBD']
        PAD_CBD_Profile_threshold = PAD_CBD_Profile[PAD_CBD_Profile['CBD_rollM'] > threshold].copy()
    
    # No data
    if len(PAD_CBD_Profile_threshold) == 0:
        return _create_null_return_ffuelmetrics(PAD_CBD_Profile)
    
    # Get number of discontinuity (FSG) of 1m or more
    shift_H = PAD_CBD_Profile_threshold['H'].shift(1)
    shift_H.iloc[0] = d / 2
    delta_layer = PAD_CBD_Profile_threshold['H'] - shift_H
    Discontinuity = delta_layer[delta_layer > 1].values
    
    # Get the FPT (Fuel Profile Type)
    Profil_Type = None
    Profil_Type_L = None
    delta_ID = None
    
    # One discontinuity = Discontinuous = Stratified = Profile 1
    if len(Discontinuity) == 1:
        delta_ID = np.where(delta_layer == Discontinuity[0])[0]
        if len(delta_ID) > 0:
            delta_ID = delta_ID[0]
        Profil_Type = 1
        Profil_Type_L = 1
    
    # If more than one discontinuities
    elif len(Discontinuity) > 1:
        # If all gaps are <= 1, keep only the first => slightly Complex
        if np.all(Discontinuity <= 1):
            Discontinuity = [Discontinuity[0]]
            delta_ID = np.where(delta_layer == Discontinuity[0])[0]
            if len(delta_ID) > 0:
                delta_ID = delta_ID[0]
            Profil_Type = 2
            Profil_Type_L = 3
        # If only one discontinuity is > 1, keep this one
        elif np.sum(Discontinuity > 1) == 1:
            Discontinuity = Discontinuity[Discontinuity > 1]
            delta_ID = np.where(delta_layer == Discontinuity[0])[0]
            if len(delta_ID) > 0:
                delta_ID = delta_ID[0]
            Profil_Type = 3
            Profil_Type_L = 3
        # If more than one discontinuity is above 1, keep the first
        elif np.sum(Discontinuity > 1) > 1:
            Discontinuity = [Discontinuity[0]]
            delta_ID = np.where(delta_layer == Discontinuity[0])[0]
            if len(delta_ID) > 0:
                delta_ID = delta_ID[0]
            Profil_Type = 4
            Profil_Type_L = 3
    
    # Get metrics
    # Profile continue = Profile 5
    if len(Discontinuity) == 0:
        CBH = 0
        FSG = 0
        Top_Fuel = np.max(PAD_CBD_Profile_threshold['H'])
        H_Bush = Top_Fuel
        continuity = 1
        Profil_Type = 5
        Profil_Type_L = 4
    
    # Profile discontinue
    if len(Discontinuity) > 0:
        # Profile discontinue without understory strata
        if np.min(PAD_CBD_Profile_threshold['H']) > 1.25:
            CBH = PAD_CBD_Profile_threshold['H'].iloc[delta_ID]
            FSG = Discontinuity[0]
            H_Bush = 0
        # Profile discontinue with understory strata
        else:
            CBH = PAD_CBD_Profile_threshold['H'].iloc[delta_ID]
            FSG = Discontinuity[0]
            H_Bush = CBH - FSG
            if Profil_Type_L == 1:
                Profil_Type_L = 2
        Top_Fuel = np.max(PAD_CBD_Profile_threshold['H'])
        continuity = 0
    
    # Get metrics (above H_PAI)
    # Need PAD column for PAI calculation
    if 'PAD' not in PAD_CBD_Profile.columns:
        # Estimate PAD from CBD if FMA is available
        if 'FMA' in PAD_CBD_Profile.columns:
            PAD_CBD_Profile['PAD'] = PAD_CBD_Profile['CBD'] / PAD_CBD_Profile['FMA']
        else:
            # Use a default FMA to estimate PAD
            default_FMA = 0.15  # Approximate value
            PAD_CBD_Profile['PAD'] = PAD_CBD_Profile['CBD'] / default_FMA
    
    profile_above = PAD_CBD_Profile[PAD_CBD_Profile['H'] >= H_PAI]
    PAI_tot = np.sum(profile_above['PAD']) * d
    
    # VCI_PAD calculation
    pad_norm = profile_above['PAD'] / np.sum(profile_above['PAD'])
    pad_norm = pad_norm[pad_norm > 0]  # Avoid log(0)
    VCI_PAD = -np.sum(pad_norm * np.log(pad_norm)) / np.log(len(pad_norm)) if len(pad_norm) > 0 else 0
    
    Height = np.max(PAD_CBD_Profile['H'])
    CBD_max = np.max(profile_above['CBD_rollM'])
    
    # Fuel load calculations
    profile_1m = PAD_CBD_Profile[(PAD_CBD_Profile['H'] > 1) & 
                                  (PAD_CBD_Profile['H'] >= CBH) & 
                                  (PAD_CBD_Profile['H'] <= Height)]
    CFL = np.sum(profile_1m['CBD_rollM']) * d
    
    profile_total = PAD_CBD_Profile[(PAD_CBD_Profile['H'] > 1) & 
                                     (PAD_CBD_Profile['H'] <= Height)]
    TFL = np.sum(profile_total['CBD_rollM']) * d
    
    if CBH == 0:
        MFL = TFL
    else:
        profile_mid = PAD_CBD_Profile[(PAD_CBD_Profile['H'] > 1) & 
                                       (PAD_CBD_Profile['H'] <= H_Bush)]
        MFL = np.sum(profile_mid['CBD_rollM']) * d
    
    FL_0_1 = np.sum(PAD_CBD_Profile[PAD_CBD_Profile['H'] <= 1]['CBD_rollM']) * d
    
    FL_1_3 = np.sum(PAD_CBD_Profile[(PAD_CBD_Profile['H'] > 1) & 
                                     (PAD_CBD_Profile['H'] <= 3)]['CBD_rollM']) * d
    
    if FSG == 0:
        GSFL = 0
    else:
        profile_gap = PAD_CBD_Profile[(PAD_CBD_Profile['H'] > H_Bush) & 
                                       (PAD_CBD_Profile['H'] <= CBH)]
        GSFL = np.sum(profile_gap['CBD_rollM']) * d
    
    # Get FMA if available
    if 'FMA' in PAD_CBD_Profile.columns:
        FMA = np.mean(PAD_CBD_Profile['FMA'])
    else:
        FMA = -1
    
    # Create metrics dictionary
    VVP_metrics = {
        'Profil_Type': Profil_Type,
        'Profil_Type_L': Profil_Type_L,
        'threshold': threshold,
        'Height': Height,
        'CBH': CBH,
        'FSG': FSG,
        'Top_Fuel': Top_Fuel,
        'H_Bush': H_Bush,
        'continuity': continuity,
        'VCI_PAD': VCI_PAD,
        'PAI_tot': PAI_tot,
        'CBD_max': CBD_max,
        'CFL': CFL,
        'TFL': TFL,
        'MFL': MFL,
        'FL_1_3': FL_1_3,
        'GSFL': GSFL,
        'FL_0_1': FL_0_1
    }
    
    if FMA != -1:
        VVP_metrics['FMA'] = FMA
    
    return (VVP_metrics, PAD_CBD_Profile)


def _create_null_return_ffuelmetrics(PAD_CBD_Profile):
    """Create a null return with -1 values for all metrics."""
    VVP_metrics = {
        'Profil_Type': -1,
        'Profil_Type_L': -1,
        'threshold': -1,
        'Height': -1,
        'CBH': -1,
        'FSG': -1,
        'Top_Fuel': -1,
        'H_Bush': -1,
        'continuity': -1,
        'VCI_PAD': -1,
        'PAI_tot': -1,
        'CBD_max': -1,
        'CFL': -1,
        'TFL': -1,
        'MFL': -1,
        'FL_1_3': -1,
        'GSFL': -1,
        'FL_0_1': -1
    }
    
    return (VVP_metrics, PAD_CBD_Profile)

