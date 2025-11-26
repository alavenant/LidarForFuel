"""
Tests for fPCpretreatment module
"""

import pytest
import tempfile
import numpy as np
from pathlib import Path
import laspy
from python.fPCpretreatment import fPCpretreatment


def create_test_las_file(x, y, z, gps_time=None, classification=None):
    """Helper function to create a test LAS file"""
    with tempfile.NamedTemporaryFile(suffix=".las", delete=False) as temp:
        las = laspy.create(point_format=3, file_version="1.4")
        las.x = x
        las.y = y
        las.z = z
        if gps_time is not None:
            las.gps_time = gps_time
        else:
            las.gps_time = np.zeros(len(x))
        # Add random classification if not provided
        if classification is not None:
            las.classification = classification
        else:
            # Random classification: 1=unclassified, 2=ground, 3-5=vegetation
            las.classification = np.random.choice([1, 2, 3, 4, 5], size=len(x))
        las.write(temp.name)
        return Path(temp.name)


def test_fPCpretreatment_basic():
    """Test basic fPCpretreatment functionality"""
    points = 100
    x = np.random.rand(points) * 1000
    y = np.random.rand(points) * 1000
    z = np.random.rand(points) * 100
    
    input_file = create_test_las_file(x, y, z)
    output_file = Path(tempfile.mktemp(suffix=".las"))
    
    try:
        las = fPCpretreatment(
            chunk=str(input_file),
        )
        
        assert las is not None, "fPCpretreatment should return a LAS object"
        assert len(las.points) > 0, "LAS should contain points"
        
        # Check that new attributes are added
        assert hasattr(las, 'LMA'), "LAS should have LMA attribute"
        assert hasattr(las, 'WD'), "LAS should have WD attribute"
        assert hasattr(las, 'Zref'), "LAS should have Zref attribute"
        
    finally:
        if input_file.exists():
            input_file.unlink()
        if output_file.exists():
            output_file.unlink()


def test_fPCpretreatment_with_custom_lma_wd():
    """Test fPCpretreatment with custom LMA and WD values"""
    points = 50
    x = np.random.rand(points) * 1000
    y = np.random.rand(points) * 1000
    z = np.random.rand(points) * 100
    
    input_file = create_test_las_file(x, y, z)
    
    try:
        H_strata_bush = 2  # Default value
        las = fPCpretreatment(
            chunk=str(input_file),
            LMA=120.6,
            WD=500,
            H_strata_bush=H_strata_bush
        )
        
        assert las is not None
        
        # For points above H_strata_bush, LMA should be the custom value
        above_bush_mask = las.z > H_strata_bush
        if np.any(above_bush_mask):
            assert np.allclose(las.LMA[above_bush_mask], 120.6), \
                "LMA should be set to custom value for points above H_strata_bush"
            assert np.allclose(las.WD[above_bush_mask], 500), \
                "WD should be set to custom value for points above H_strata_bush"
        
        # For points at or below H_strata_bush, LMA should be LMA_bush (default 140)
        below_bush_mask = las.z <= H_strata_bush
        if np.any(below_bush_mask):
            assert np.allclose(las.LMA[below_bush_mask], 140), \
                "LMA should be LMA_bush (default 140) for points at or below H_strata_bush"
            assert np.allclose(las.WD[below_bush_mask], 591), \
                "WD should be WD_bush (default 591) for points at or below H_strata_bush"
        
    finally:
        if input_file.exists():
            input_file.unlink()


def test_fPCpretreatment_season_filter():
    """Test fPCpretreatment with season filter"""
    points = 100
    x = np.random.rand(points) * 1000
    y = np.random.rand(points) * 1000
    z = np.random.rand(points) * 100
    # Create gps_time that corresponds to summer months (May-October)
    gps_time = np.linspace(0, 365*24*3600, points)  # One year of data
    
    input_file = create_test_las_file(x, y, z, gps_time)
    
    try:
        las = fPCpretreatment(
            chunk=str(input_file),
            season_filter=[5, 6, 7, 8, 9, 10],  # May to October
            start_date="2023-01-01 00:00:00"
        )
        
        # Should have filtered points (exact number depends on dates)
        assert las is not None
        
    finally:
        if input_file.exists():
            input_file.unlink()


def test_fPCpretreatment_height_filter():
    """Test fPCpretreatment with height filter"""
    points = 100
    x = np.random.rand(points) * 1000
    y = np.random.rand(points) * 1000
    z = np.random.rand(points) * 200  # Some points above 60m
    
    input_file = create_test_las_file(x, y, z)
    
    try:
        las = fPCpretreatment(
            chunk=str(input_file),
            Height_filter=60
        )
        
        assert las is not None
        if len(las.points) > 0:
            assert np.all(las.z < 60), "All points should be below height filter"
        
    finally:
        if input_file.exists():
            input_file.unlink()


if __name__ == "__main__":
    pytest.main([__file__])

