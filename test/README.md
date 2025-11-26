# Tests for LidarForFuel Python Package

This directory contains tests for the Python implementation of LidarForFuel.

## Running Tests

### Install test dependencies

```bash
pip install pytest pytest-cov
```

Or using conda:

```bash
conda install pytest pytest-cov
```

### Run all tests

```bash
pytest test/
```

### Run specific test file

```bash
pytest test/test_fPCpretreatment.py
pytest test/test_fCBDprofile_fuelmetrics.py
pytest test/test_ffuelmetrics.py
```

### Run with coverage

```bash
pytest test/ --cov=python --cov-report=html
```

### Run with verbose output

```bash
pytest test/ -v
```

## Test Structure

- `test_fPCpretreatment.py`: Tests for point cloud preprocessing
- `test_fCBDprofile_fuelmetrics.py`: Tests for fuel metrics computation from point clouds
- `test_ffuelmetrics.py`: Tests for fuel metrics computation from profiles
- `conftest.py`: Shared fixtures and pytest configuration

## Writing New Tests

When adding new tests:

1. Follow the naming convention: `test_*.py` for test files, `test_*` for test functions
2. Use fixtures from `conftest.py` when possible
3. Clean up temporary files in `finally` blocks or using fixtures
4. Add docstrings to test functions explaining what they test

## Requirements

Tests require:
- pytest
- numpy
- pandas
- laspy (for point cloud tests)

All dependencies are listed in `../python/requirements.txt`.

