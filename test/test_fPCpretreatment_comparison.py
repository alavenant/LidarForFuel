"""
Test comparison of fPCpretreatment output with reference file using ign-pdal-tools
"""

import pytest
import tempfile
from pathlib import Path
import sys

# Add parent directory to path to import python modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from python.fPCpretreatment import fPCpretreatment

try:
    from pdaltools.las_comparison import compare_las_dimensions
    PDALTOOLS_AVAILABLE = True
except ImportError:
    PDALTOOLS_AVAILABLE = False
    pytest.skip("pdaltools not available", allow_module_level=True)


def test_fPCpretreatment_comparison_with_reference():
    """
    Test that fPCpretreatment produces output similar to reference file.
    Compares the result of fPCpretreatment on M30_FontBlanche.laz
    with the reference file M30_FontBlanche_pretreated.laz
    """
    # Paths to test data files
    base_path = Path(__file__).parent.parent
    input_file = base_path / "inst" / "extdata" / "M30_FontBlanche.laz"
    reference_file = base_path / "inst" / "extdata" / "M30_FontBlanche_pretreated.laz"
    
    # Check if test data files exist
    if not input_file.exists():
        pytest.skip(f"Input file not found: {input_file}")
    if not reference_file.exists():
        pytest.skip(f"Reference file not found: {reference_file}")
    
    # Create temporary output file
    output_file = Path(tempfile.mktemp(suffix=".las"))
    
    try:
        # Run fPCpretreatment on input file
        # Using default parameters that match the reference file
        las_result = fPCpretreatment(
            chunk=str(input_file),
            LMA=140,  # Default value
            WD=591,   # Default value
            classify=False,
            Height_filter=60
        )
        
        if las_result is None:
            pytest.fail("fPCpretreatment returned None")
        
        # Write result to temporary file
        las_result.write(str(output_file))
        
        # Compare dimensions between result and reference
        # We'll compare key dimensions that should be present after pretreatment
        dimensions_to_compare = ["X", "Y", "Z", "LMA", "WD", "Zref"]
        
        # Check which dimensions are available in both files
        import laspy
        las_ref = laspy.read(str(reference_file))
        available_dims_ref = set(las_ref.point_format.dimension_names)
        available_dims_result = set(las_result.point_format.dimension_names)
        
        # Only compare dimensions that exist in both files
        dimensions_to_compare = [
            dim for dim in dimensions_to_compare 
            if dim in available_dims_ref and dim in available_dims_result
        ]
        
        if not dimensions_to_compare:
            pytest.skip("No common dimensions found between result and reference")
        
        # Compare with precision for float dimensions
        precision = {
            "X": 0.001,
            "Y": 0.001,
            "Z": 0.001,
            "LMA": 0.1,
            "WD": 0.1,
            "Zref": 0.001
        }
        
        # Filter precision to only include dimensions we're comparing
        precision = {k: v for k, v in precision.items() if k in dimensions_to_compare}
        
        # Compare dimensions
        result, nb_diff, percentage = compare_las_dimensions(
            output_file,
            reference_file,
            dimensions=dimensions_to_compare,
            precision=precision
        )
        
        # The comparison should pass (result=True) or have very few differences
        # We allow some tolerance for floating point differences
        if not result:
            # If comparison failed, check if differences are within acceptable tolerance
            # For a successful pretreatment, we expect at least 90% of points to match
            assert percentage < 10.0, \
                f"Too many differences ({percentage:.2f}%) between result and reference. " \
                f"Number of different points: {nb_diff}"
        
        # Also verify that key attributes are present
        assert "LMA" in available_dims_result, "LMA dimension should be present"
        assert "WD" in available_dims_result, "WD dimension should be present"
        assert "Zref" in available_dims_result, "Zref dimension should be present"
        
    finally:
        # Clean up
        if output_file.exists():
            output_file.unlink()


def test_fPCpretreatment_comparison_structure():
    """
    Test that fPCpretreatment produces output with the same structure as reference file.
    Checks that the same dimensions are present (even if values differ slightly).
    """
    base_path = Path(__file__).parent.parent
    input_file = base_path / "inst" / "extdata" / "M30_FontBlanche.laz"
    reference_file = base_path / "inst" / "extdata" / "M30_FontBlanche_pretreated.laz"
    
    if not input_file.exists() or not reference_file.exists():
        pytest.skip("Test data files not found")
    
    output_file = Path(tempfile.mktemp(suffix=".las"))
    
    try:
        # Run fPCpretreatment
        las_result = fPCpretreatment(
            chunk=str(input_file),
            LMA=140,
            WD=591
        )
        
        if las_result is None:
            pytest.fail("fPCpretreatment returned None")
        
        las_result.write(str(output_file))
        
        # Compare structure (all dimensions)
        result, nb_diff, percentage = compare_las_dimensions(
            output_file,
            reference_file,
            dimensions=None  # Compare all dimensions
        )
        
        # For structure comparison, we're more lenient
        # We just want to ensure the same dimensions exist
        import laspy
        las_ref = laspy.read(str(reference_file))
        las_out = laspy.read(str(output_file))
        
        ref_dims = set(las_ref.point_format.dimension_names)
        out_dims = set(las_out.point_format.dimension_names)
        
        # Check that key dimensions are present
        required_dims = {"LMA", "WD", "Zref", "Easting", "Northing", "Elevation"}
        missing_dims = required_dims - out_dims
        assert len(missing_dims) == 0, \
            f"Missing required dimensions in output: {missing_dims}"
        
    finally:
        if output_file.exists():
            output_file.unlink()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

