"""
Fuel metrics LiDAR

This module provides functions to compute PAD and CBD profiles from pretreated ALS point clouds
and obtain fuel metrics from them.
"""

import numpy as np
import pandas as pd
import warnings
from typing import Union, List, Tuple, Optional
try:
    import laspy
except ImportError:
    laspy = None


def fCBDprofile_fuelmetrics(
    datatype="Pixel",
    X=None, Y=None, Z=None, Zref=None,
    ReturnNumber=None,
    Easting=None, Northing=None, Elevation=None,
    LMA=None,
    gpstime=None,
    Height_Cover=2,
    threshold=0.02,
    scanning_angle=True,
    use_cover=False,
    WD=None,
    limit_N_points=400,
    limit_flightheight=800,
    limit_vegetationheight=0.1,
    H_PAI=0,
    omega=0.77,
    d=1,
    G=0.5
):
    """
    Compute PAD and CBD profiles from a pretreated ALS point cloud and obtain fuel metrics.
    
    Parameters
    ----------
    datatype : str or laspy.LasData, default "Pixel"
        Either "Pixel" or directly a las/laz file. "Pixel" if used with pixel_metric function
        to map fuel metrics. Or a laspy.LasData object if a plot point cloud only needs to be computed.
    X, Y, Z, Zref : array-like, optional
        Coordinates of a point cloud (Z being the normalized Z coordinate and Zref the original one)
    ReturnNumber : array-like, optional
        The return number of a pulse. Used for calculating cover.
    Easting, Northing, Elevation : array-like, optional
        Coordinates of the plane associated to each point
    LMA : array-like or numeric, optional
        Leaf mass area in g/m² associated to each point or a generic value
    gpstime : array-like, optional
        GPS time of the point cloud. Only used to retrieve date of scanning
    Height_Cover : numeric, default 2
        The height from which the canopy cover should be estimated.
    threshold : numeric or str, default 0.02
        Bulk density critical threshold used to discriminate different strata limits.
        Either numeric: a bulk density value (in kg/m³) or str: a percentage of maximum CBD value.
    scanning_angle : bool, default True
        Use the scanning angle computed from the trajectories to estimate cos(theta).
        If false: cos(theta) = 1
    use_cover : bool, default False
        Use cover for PAD estimates
    WD : array-like or numeric, optional
        Wood density associated to each point or a generic value
    limit_N_points : numeric, default 400
        Minimum number of points in the pixel/plot for computing profiles & metrics.
    limit_flightheight : numeric, default 800
        Flight height above canopy in m. If the flight height is lower than limit_flightheight,
        bulk density profile is not computed.
    limit_vegetationheight : numeric, default 0.1
        Vegetation height below which bulk density profile should not be computed.
    H_PAI : numeric, default 0
        Height from which PAI, VCI and entropy is estimated. Default is 0 (means from ground to top).
    omega : numeric, default 0.77
        Clumping factor. Default is 0.77. One means "no clumping".
    d : numeric, default 1
        Depth of the strata in meter to compute the profile
    G : numeric, default 0.5
        Leaf projection ratio.
        
    Returns
    -------
    If datatype = "Pixel": dict
        Dictionary with 173 keys corresponding to metrics and bulk density profile values
        per strata of depth d.
    If datatype is a las: tuple
        A tuple of two elements:
        1) A dict with all fuel metrics
        2) A pandas DataFrame with the PAD and CBD profile values (columns: H, PAD, CBD)
    """
    # Extract data from LAS object if provided
    if laspy is not None and isinstance(datatype, laspy.LasData):
        X = datatype.x
        Y = datatype.y
        Z = datatype.z
        if hasattr(datatype, 'Zref'):
            Zref = datatype.Zref
        else:
            Zref = Z.copy()
        ReturnNumber = datatype.return_number
        if hasattr(datatype, 'Easting'):
            Easting = datatype.Easting
        if hasattr(datatype, 'Northing'):
            Northing = datatype.Northing
        if hasattr(datatype, 'Elevation'):
            Elevation = datatype.Elevation
        if hasattr(datatype, 'LMA'):
            LMA = datatype.LMA
        if hasattr(datatype, 'WD'):
            WD = datatype.WD
        if hasattr(datatype, 'gps_time'):
            gpstime = datatype.gps_time
        elif hasattr(datatype, 'Time'):
            gpstime = datatype.Time
    
    # Convert to numpy arrays
    X = np.asarray(X)
    Y = np.asarray(Y)
    Z = np.asarray(Z)
    if Zref is None:
        Zref = Z.copy()
    else:
        Zref = np.asarray(Zref)
    ReturnNumber = np.asarray(ReturnNumber)
    Easting = np.asarray(Easting) if Easting is not None else X.copy()
    Northing = np.asarray(Northing) if Northing is not None else Y.copy()
    Elevation = np.asarray(Elevation) if Elevation is not None else Z.copy()
    LMA = np.asarray(LMA) if LMA is not None else np.full(len(Z), 140.0)
    WD = np.asarray(WD) if WD is not None else np.full(len(Z), 591.0)
    if gpstime is not None:
        gpstime = np.asarray(gpstime)
        date = np.mean(gpstime)
    else:
        date = 0
    
    # Check minimum number of points
    if len(Z) < limit_N_points:
        warnings.warn(
            "NULL return: The number of points < limit_N_points: "
            "check your tile or your point cloud"
        )
        return _create_null_return(date, datatype)
    
    # Get cover
    cover = np.sum((ReturnNumber[Z > Height_Cover] == 1)) / np.sum(ReturnNumber == 1) if np.sum(ReturnNumber == 1) > 0 else 0
    cover_4 = np.sum((ReturnNumber[Z > 4] == 1)) / np.sum(ReturnNumber == 1) if np.sum(ReturnNumber == 1) > 0 else 0
    cover_6 = np.sum((ReturnNumber[Z > 6] == 1)) / np.sum(ReturnNumber == 1) if np.sum(ReturnNumber == 1) > 0 else 0
    
    # Get PAD and CBD profile
    # Create a sequence to make strata
    seq_layer = np.concatenate([[np.min(Z)], np.arange(0, np.max(Z) + d, d), [np.max(Z)]])
    
    # Histogram to get number of returns in strata
    Ni, _ = np.histogram(Z, bins=seq_layer)
    N = np.cumsum(Ni)
    NRD = Ni / N
    NRD = np.nan_to_num(NRD, nan=0.0)
    
    # NRD estimation (equations 23 and 24 of Pimont et al 2018)
    NRD[NRD == 0] = (Ni[NRD == 0] + 1) / (N[NRD == 0] + 2)
    NRD[NRD == 1] = (Ni[NRD == 1] + 1) / (N[NRD == 1] + 2)
    
    # Gap fraction estimation
    Gf = 1 - NRD
    
    # Calculate component of vector U (plane -> point) to take into account scanning angle
    if scanning_angle:
        norm_U = np.sqrt((X - Easting)**2 + (Y - Northing)**2 + (Zref - Elevation)**2)
        Nx_U = np.abs((X - Easting) / norm_U)
        Ny_U = np.abs((Y - Northing) / norm_U)
        Nz_U = np.abs((Zref - Elevation) / norm_U)
    else:
        Nz_U = np.ones_like(Z)
        norm_U = np.full_like(Z, 999999)
    
    # Exception if mean of norm_U < limit_flightheight
    if np.mean(norm_U[~np.isnan(norm_U)]) < limit_flightheight:
        warnings.warn(
            "NULL return: limit_flightheight below the threshold. "
            "Check your trajectory and avoid using scanning_angle mode if the trajectory is uncertain"
        )
        return _create_null_return(date, datatype)
    
    # Remove the bottom & top value of seq and add d/2 to get the middle height of strata
    seq_layer = seq_layer[1:-1] + (d / 2)
    
    # cos theta takes into account scanning angle
    cos_theta = np.mean(np.abs(Nz_U[~np.isnan(Nz_U)]))
    
    # Plant area density calculation (FAD -> fuel area density: leaves + twigs)
    PAD = -(np.log(Gf) * cos_theta / (G * omega) / d)
    
    if use_cover:
        PAD = (-np.log(1 - Ni / (N * cover)) / (G * omega * (d / cos_theta))) * cover
        if Height_Cover >= np.max(seq_layer):
            PAD = -(np.log(Gf) * cos_theta / (G * omega) / d)
            warnings.warn("Cover method was not used as Height_Cover > Vegetation Height")
        if cover == 0:
            PAD = -(np.log(Gf) * cos_theta / (G * omega) / d)
            warnings.warn("Cover method was not used as Cover = 0")
    
    # Standard deviation of PAD
    SD_PAD = (2 / d) * np.sqrt(NRD / (N * (1 - NRD)))
    SD_PAD[NRD == 1] = (2 / d) * np.sqrt(2 + 1 / N[NRD == 1])
    
    # LMA from g/cm² to kg/m²
    LMA_mean = np.mean(LMA[~np.isnan(LMA)]) / 1000
    if np.isnan(LMA_mean):
        LMA_mean = 0.15
    
    # Partition of fuel surface (fine branch vs leaves)
    # Wood density (kg/m³)
    WD_mean = np.mean(WD[~np.isnan(WD)])
    
    # Surface volume ratio (SVR: m²/m³) for 4mm diameter twigs
    SVR = 1 / 0.002
    # Wood mass area (WMA)
    WMA = WD_mean / SVR
    
    # Partition of wood and leaves (M. Soma phd thesis data)
    partW = 0.51
    partL = 0.49
    
    # Fuel mass area
    FMA = 1 / ((partW / WMA) + (partL / LMA_mean))
    
    # CBD in kg/m³
    CBD = PAD * FMA
    
    # Table with strata height, PAD and CBD
    PAD_CBD_Profile = pd.DataFrame({
        'H': seq_layer,
        'PAD': PAD[1:],
        'CBD': CBD[1:],
        'SD_PAD': SD_PAD[1:],
        'NRD': NRD[1:],
        'Ni': Ni[1:],
        'N': N[1:]
    })
    
    # Work on profile to get FPT and fuel metrics
    # Define threshold when threshold is a proportion of CBD max
    if isinstance(threshold, str) and '%' in threshold:
        threshold_prop = float(threshold.replace('%', '')) / 100
        Height_CBD = max(np.max(PAD_CBD_Profile['H']) / 3, 2)
        threshold = np.max(PAD_CBD_Profile[PAD_CBD_Profile['H'] > Height_CBD]['CBD']) * threshold_prop
    else:
        threshold = float(threshold)
    
    # No data above limit_vegetationheight
    if np.max(PAD_CBD_Profile['H']) < limit_vegetationheight:
        warnings.warn(f"NULL (-1) return: no data above {limit_vegetationheight}m height")
        return _create_null_return(date, datatype)
    
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
        return _create_null_return(date, datatype)
    
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
    profile_above = PAD_CBD_Profile[PAD_CBD_Profile['H'] >= H_PAI]
    PAI_tot = np.sum(profile_above['PAD']) * d
    
    # VCI_PAD calculation
    pad_norm = profile_above['PAD'] / np.sum(profile_above['PAD'])
    pad_norm = pad_norm[pad_norm > 0]  # Avoid log(0)
    VCI_PAD = -np.sum(pad_norm * np.log(pad_norm)) / np.log(len(pad_norm)) if len(pad_norm) > 0 else 0
    
    # VCI_lidr and entropy_lidr (simplified versions)
    Z_above = Z[Z >= H_PAI]
    if len(Z_above) > 0:
        z_max = np.max(Z)
        # Simplified VCI calculation
        VCI_lidr = np.std(Z_above) / z_max if z_max > 0 else 0
        # Simplified entropy calculation
        hist, _ = np.histogram(Z_above, bins=10)
        hist_norm = hist / np.sum(hist)
        hist_norm = hist_norm[hist_norm > 0]
        entropy_lidr = -np.sum(hist_norm * np.log(hist_norm)) if len(hist_norm) > 0 else 0
    else:
        VCI_lidr = 0
        entropy_lidr = 0
    
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
        'VCI_lidr': VCI_lidr,
        'entropy_lidr': entropy_lidr,
        'PAI_tot': PAI_tot,
        'CBD_max': CBD_max,
        'CFL': CFL,
        'TFL': TFL,
        'MFL': MFL,
        'FL_1_3': FL_1_3,
        'GSFL': GSFL,
        'FL_0_1': FL_0_1,
        'FMA': FMA,
        'date': date,
        'Cover': cover,
        'Cover_4': cover_4,
        'Cover_6': cover_6
    }
    
    # Add CBD profile values (up to 150 layers)
    VVP_metrics_CBD = np.full(150, -1.0)
    cbd_values = PAD_CBD_Profile['CBD_rollM'].values
    VVP_metrics_CBD[:len(cbd_values)] = cbd_values
    
    for i in range(150):
        VVP_metrics[f'CBD_{i+1}'] = VVP_metrics_CBD[i]
    
    # Return format depends on datatype
    if laspy is not None and isinstance(datatype, laspy.LasData):
        return (VVP_metrics, PAD_CBD_Profile)
    else:
        return VVP_metrics


def _create_null_return(date, datatype):
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
        'VCI_lidr': -1,
        'entropy_lidr': -1,
        'PAI_tot': -1,
        'CBD_max': -1,
        'CFL': -1,
        'TFL': -1,
        'MFL': -1,
        'FL_1_3': -1,
        'GSFL': -1,
        'FL_0_1': -1,
        'FMA': -1,
        'date': date,
        'Cover': -1,
        'Cover_4': -1,
        'Cover_6': -1
    }
    
    for i in range(150):
        VVP_metrics[f'CBD_{i+1}'] = -1.0
    
    PAD_CBD_Profile = None
    
    if laspy is not None and isinstance(datatype, laspy.LasData):
        return (VVP_metrics, PAD_CBD_Profile)
    else:
        return VVP_metrics

