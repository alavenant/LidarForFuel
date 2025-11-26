"""
Pytest configuration and shared fixtures
"""

import pytest
import numpy as np
import pandas as pd
import tempfile
from pathlib import Path
import laspy


@pytest.fixture
def sample_point_cloud():
    """Create a sample point cloud for testing"""
    points = 100
    x = np.random.rand(points) * 1000
    y = np.random.rand(points) * 1000
    z = np.random.rand(points) * 50
    
    return x, y, z


@pytest.fixture
def temp_las_file(sample_point_cloud):
    """Create a temporary LAS file"""
    x, y, z = sample_point_cloud
    
    with tempfile.NamedTemporaryFile(suffix=".las", delete=False) as temp:
        las = laspy.create(point_format=3, file_version="1.4")
        las.x = x
        las.y = y
        las.z = z
        las.gps_time = np.zeros(len(x))
        las.write(temp.name)
        file_path = Path(temp.name)
    
    yield file_path
    
    # Cleanup
    if file_path.exists():
        file_path.unlink()


@pytest.fixture
def sample_profile():
    """Create a sample profile DataFrame for testing"""
    return pd.DataFrame({
        'H': np.arange(0, 10, 0.5),
        'BD': np.random.rand(20) * 0.1 + 0.01,
        'PAD': np.random.rand(20) * 2.0 + 0.1
    })

