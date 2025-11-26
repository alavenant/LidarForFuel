"""
LidarForFuel - Python Package

This package provides tools to compute fuel metrics from LiDAR point clouds.
The method is fully described in the paper Martin-Ducup et al 2024.

The package contains three main modules:
- fPCpretreatment: Preprocessing of LAS/LAZ files
- fCBDprofile_fuelmetrics: Compute fuel metrics from pretreated point clouds
- ffuelmetrics: Compute fuel metrics from bulk density profiles
"""

from .fPCpretreatment import fPCpretreatment
from .fCBDprofile_fuelmetrics import fCBDprofile_fuelmetrics
from .ffuelmetrics import ffuelmetrics

__version__ = "1.0.0"
__all__ = ['fPCpretreatment', 'fCBDprofile_fuelmetrics', 'ffuelmetrics']

