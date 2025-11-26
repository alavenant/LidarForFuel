"""
Point cloud pre-treatment for using fCBDprofile_fuelmetrics in pixels

This module provides functions for preprocessing LAS/LAZ files for fuel metrics computation.
"""

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from pathlib import Path
import warnings
import argparse
import sys
try:
    import laspy
except ImportError:
    laspy = None
try:
    import rasterio
    from rasterio.mask import mask
except ImportError:
    rasterio = None
try:
    from scipy.spatial import cKDTree
except ImportError:
    cKDTree = None


def fPCpretreatment(chunk, classify=False, LMA=140, WD=591, WD_bush=591, 
                    LMA_bush=140, H_strata_bush=2, Height_filter=60,
                    start_date="2011-09-14 00:00:00", season_filter=None,
                    deviation_days=None, plot_hist_days=False):
    """
    Preprocess LAS/LAZ files for use in fCBDprofile_fuelmetrics.
    
    Parameters
    ----------
    chunk : str or laspy.LasData
        Path to a LAS/LAZ file or a laspy.LasData object
    classify : bool, default False
        Make a ground classification. Only if the original point cloud is not classified
    LMA : str or numeric, default 140
        Path to a LMA map (.tif) or a single LMA value in g/m²
    WD : str or numeric, default 591
        Path to a WD map (.tif) or a single WD value in kg/m³
    LMA_bush : str or numeric, default 140
        Similar to LMA but for the understorey strata 0 to 2m
    WD_bush : str or numeric, default 591
        Similar to WD but for the understorey strata 0 to 2m
    H_strata_bush : numeric, default 2
        Height of the strata to consider for separating LMA and WD between canopy and bush
    Height_filter : numeric, default 60
        Height limit to remove noise points
    start_date : str, default "2011-09-14 00:00:00"
        The absolute starting date to retrieve date from relative gpstime
    season_filter : list of int, default None (all months)
        A list of integers for months to keep (e.g., [5,6,7,8,9,10] for May to October)
    deviation_days : numeric, default None
        Maximum number of days tolerated between the acquisition in a given point cloud
    plot_hist_days : bool, default False
        Should the histogram of dates of acquisition be displayed
        
    Returns
    -------
    laspy.LasData or None
        Normalized point cloud with new attributes needed to run fCBDprofile_fuelmetrics
    """
    if laspy is None:
        raise ImportError("laspy is required for this function. Install it with: pip install laspy")
    
    # Read chunk with proper error handling for LAZ files
    if isinstance(chunk, str):
        try:
            las = laspy.read(chunk)
        except laspy.errors.LaspyException as e:
            error_msg = str(e)
            if "No LazBackend" in error_msg or "cannot decompress" in error_msg or "LazBackend" in error_msg:
                raise ImportError(
                    f"LAZ decompression error: {error_msg}\n\n"
                    "LAZ files require a decompression backend. Please install one of the following:\n"
                    "  - lazrs (recommended, fastest): pip install lazrs\n"
                    "  - laszip (alternative): pip install laszip\n"
                    "Or install laspy with a backend: pip install 'laspy[lazrs]'\n\n"
                    "After installation, try running the script again."
                ) from e
            else:
                # Re-raise other laspy errors as-is
                raise
        except Exception as e:
            # Handle other potential errors (file not found, etc.)
            raise
    else:
        las = chunk
    
    start_date = datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S")
    
    # Convert gpstime to date
    if hasattr(las, 'gps_time'):
        gpstime = las.gps_time
    else:
        # If no gps_time, create a dummy array
        gpstime = np.zeros(len(las.points))
    
    new_date = np.array([start_date + timedelta(seconds=float(t)) for t in gpstime])
    
    # Test season
    if season_filter is None:
        season_filter = list(range(1, 13))  # All months
    
    months_acquisition = np.array([d.month for d in new_date])
    pts_summer = np.isin(months_acquisition, season_filter)
    
    proportions_of_winter_point = (1 - np.sum(pts_summer) / len(months_acquisition)) * 100
    
    # Filter points by season
    las.points = las.points[pts_summer]
    new_date = new_date[pts_summer]
    months_acquisition = months_acquisition[pts_summer]
    
    if not np.all(np.isin(months_acquisition, season_filter)):
        if plot_hist_days:
            import matplotlib.pyplot as plt
            plt.hist(new_date, bins='auto')
            plt.title("Histogram of acquisition date")
            plt.xlabel("Date of acquisition")
            plt.show()
        month_names = [datetime(2000, m, 1).strftime("%b") for m in season_filter]
        warnings.warn(
            f"Careful {proportions_of_winter_point:.1f}% of the returns were excluded "
            f"because they were sampled outside of the chosen season "
            f"(Month: {' '.join(month_names)})"
        )
    
    if len(las.points) == 0:
        print("No points in the point cloud")
        return None
    
    # Deviation days filter
    if deviation_days is not None and np.isfinite(deviation_days):
        # Create histogram of dates
        dates_only = np.array([d.date() for d in new_date])
        unique_dates, counts = np.unique(dates_only, return_counts=True)
        
        if len(unique_dates) > deviation_days:
            max_idx = np.argmax(counts)
            max_count_date = unique_dates[max_idx]
            
            # Find dates within deviation_days
            date_diffs = np.abs([(d - max_count_date).days for d in unique_dates])
            good_date_indices = np.where(date_diffs <= deviation_days)[0]
            good_dates = unique_dates[good_date_indices]
            
            # Filter points
            keep_mask = np.isin(dates_only, good_dates)
            las.points = las.points[keep_mask]
            new_date = new_date[keep_mask]
            
            percentage_point_remove = (1 - len(las.points) / len(dates_only)) * 100
            warnings.warn(
                f"Careful {percentage_point_remove:.1f}% of the returns were removed "
                f"because they had a deviation of days around the most abundant date "
                f"greater than your threshold ({deviation_days} days)."
            )
    
    # Track sensor trajectory (simplified version)
    # In R, this uses lidR::track_sensor which is complex
    # For now, we'll use a simplified approach
    try:
        # Try to get ground points for trajectory estimation
        if hasattr(las, 'classification'):
            ground_mask = las.classification == 2  # Ground class
        else:
            ground_mask = np.ones(len(las.points), dtype=bool)
        
        if np.sum(ground_mask) == 0:
            warnings.warn("Only ground points in the tile. NULL returned")
            print("Only ground points in the tile. NULL returned")
            return None
        
        # Simplified trajectory: use mean coordinates of ground points
        ground_points = las.points[ground_mask]
        traj = pd.DataFrame({
            'Easting': [np.mean(las.x[ground_mask])],
            'Northing': [np.mean(las.y[ground_mask])],
            'Elevation': [np.mean(las.z[ground_mask]) + 1400],  # Approximate flight height
            'Time': [np.mean(gpstime[ground_mask])]
        })
        
    except Exception as e:
        warnings.warn(f"Error in trajectory estimation: {e}. Using simplified approach.")
        traj = pd.DataFrame({
            'Easting': [np.mean(las.x)],
            'Northing': [np.mean(las.y)],
            'Elevation': [np.mean(las.z) + 1400],
            'Time': [np.mean(gpstime)]
        })
    
    # Find closest gpstime between traj and las
    if cKDTree is not None:
        tree = cKDTree(traj[['Time']].values)
        _, nn_indices = tree.query(gpstime.reshape(-1, 1), k=1)
    else:
        # Fallback: use nearest neighbor search
        nn_indices = np.array([
            np.argmin(np.abs(traj['Time'].values - t)) 
            for t in gpstime
        ])
    
    # Add trajectory attributes to points
    las.add_extra_dim(laspy.ExtraBytesParams(name="Easting", type=np.float64))
    las.add_extra_dim(laspy.ExtraBytesParams(name="Northing", type=np.float64))
    las.add_extra_dim(laspy.ExtraBytesParams(name="Elevation", type=np.float64))
    las.add_extra_dim(laspy.ExtraBytesParams(name="Time", type=np.float64))
    
    las.Easting = traj.iloc[nn_indices]['Easting'].values
    las.Northing = traj.iloc[nn_indices]['Northing'].values
    las.Elevation = traj.iloc[nn_indices]['Elevation'].values
    las.Time = traj.iloc[nn_indices]['Time'].values
    
    # Ground classification (if requested)
    if classify:
        # This would require a ground classification algorithm
        # For now, we'll skip it as it's complex
        warnings.warn("Ground classification requested but not implemented. "
                     "Please classify your point cloud beforehand.")
    
    # LMA and WD assignment
    if isinstance(LMA, (int, float)):
        las.add_extra_dim(laspy.ExtraBytesParams(name="LMA", type=np.float64))
        las.LMA = np.full(len(las.points), LMA)
    elif isinstance(LMA, str) and rasterio is not None:
        # Load LMA map and extract values
        with rasterio.open(LMA) as src:
            coords = list(zip(las.x, las.y))
            values = [v[0] for v in src.sample(coords)]
        las.add_extra_dim(laspy.ExtraBytesParams(name="LMA", type=np.float64))
        las.LMA = np.array(values)
    else:
        las.add_extra_dim(laspy.ExtraBytesParams(name="LMA", type=np.float64))
        las.LMA = np.full(len(las.points), 140.0)
    
    if isinstance(WD, (int, float)):
        las.add_extra_dim(laspy.ExtraBytesParams(name="WD", type=np.float64))
        las.WD = np.full(len(las.points), WD)
    elif isinstance(WD, str) and rasterio is not None:
        # Load WD map and extract values
        with rasterio.open(WD) as src:
            coords = list(zip(las.x, las.y))
            values = [v[0] for v in src.sample(coords)]
        las.add_extra_dim(laspy.ExtraBytesParams(name="WD", type=np.float64))
        las.WD = np.array(values)
    else:
        las.add_extra_dim(laspy.ExtraBytesParams(name="WD", type=np.float64))
        las.WD = np.full(len(las.points), 591.0)
    
    # Normalize height (simplified - would need DTM generation)
    # For now, we assume Z is already normalized or we use a simple approach
    if hasattr(las, 'z'):
        # Store original Z
        las.add_extra_dim(laspy.ExtraBytesParams(name="Zref", type=np.float64))
        las.Zref = las.z.copy()
        
        # Simple normalization: subtract minimum Z of ground points
        if hasattr(las, 'classification'):
            ground_z = las.z[las.classification == 2]
            if len(ground_z) > 0:
                las.z = las.z - np.min(ground_z)
        else:
            # Use lowest 2% of points as ground approximation
            z_sorted = np.sort(las.z)
            ground_level = z_sorted[int(len(z_sorted) * 0.02)]
            las.z = las.z - ground_level
    
    # Filter points
    if hasattr(las, 'classification'):
        keep_mask = (las.classification <= 5) & (las.z < Height_filter)
    else:
        keep_mask = las.z < Height_filter
    
    las.points = las.points[keep_mask]
    
    # Apply bush LMA and WD for low strata
    if hasattr(las, 'z'):
        bush_mask = las.z <= H_strata_bush
        if np.any(bush_mask):
            las.LMA[bush_mask] = LMA_bush
            las.WD[bush_mask] = WD_bush
    
    return las


def main():
    """
    Main function for command-line interface.
    """
    parser = argparse.ArgumentParser(
        description='Preprocess LAS/LAZ files for fuel metrics computation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default parameters
  python fPCpretreatment.py -i input.laz -o output.laz

  # With custom LMA and WD values
  python fPCpretreatment.py -i input.laz -o output.laz --LMA 120.6 --WD 500

  # With LMA and WD maps
  python fPCpretreatment.py -i input.laz -o output.laz --LMA lma_map.tif --WD wd_map.tif

  # Filter by season (May to October)
  python fPCpretreatment.py -i input.laz -o output.laz --season-filter 5 6 7 8 9 10

  # With deviation days filter
  python fPCpretreatment.py -i input.laz -o output.laz --deviation-days 5
        """
    )
    
    parser.add_argument('-i', '--input', type=str, required=True,
                       help='Path to input LAS/LAZ file')
    parser.add_argument('-o', '--output', type=str, required=True,
                       help='Path to output LAS/LAZ file')
    parser.add_argument('--classify', action='store_true',
                       help='Make a ground classification (default: False)')
    parser.add_argument('--LMA', type=str, default='140',
                       help='LMA value (numeric) or path to LMA map (.tif). Default: 140')
    parser.add_argument('--WD', type=str, default='591',
                       help='WD value (numeric) or path to WD map (.tif). Default: 591')
    parser.add_argument('--LMA-bush', type=float, default=140,
                       help='LMA for understorey strata 0 to 2m. Default: 140')
    parser.add_argument('--WD-bush', type=float, default=591,
                       help='WD for understorey strata 0 to 2m. Default: 591')
    parser.add_argument('--H-strata-bush', type=float, default=2,
                       help='Height of strata to consider for separating LMA/WD. Default: 2')
    parser.add_argument('--Height-filter', type=float, default=60,
                       help='Height limit to remove noise points. Default: 60')
    parser.add_argument('--start-date', type=str, default='2011-09-14 00:00:00',
                       help='Starting date for GPS time conversion. Default: "2011-09-14 00:00:00"')
    parser.add_argument('--season-filter', type=int, nargs='+', default=None,
                       help='Months to keep (e.g., 5 6 7 8 9 10 for May-October). Default: all months')
    parser.add_argument('--deviation-days', type=float, default=None,
                       help='Maximum days deviation from most abundant date. Default: None (no filter)')
    parser.add_argument('--plot-hist-days', action='store_true',
                       help='Display histogram of acquisition dates')
    
    args = parser.parse_args()
    
    # Get input and output files from arguments
    input_file = args.input
    output_file = args.output
    
    # Convert LMA and WD to appropriate types
    try:
        LMA = float(args.LMA)
    except ValueError:
        LMA = args.LMA  # Assume it's a file path
    
    try:
        WD = float(args.WD)
    except ValueError:
        WD = args.WD  # Assume it's a file path
    
    # Check if input file exists
    input_path = Path(input_file)
    if not input_path.exists():
        print(f"Error: Input file '{input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Processing: {input_file}")
    print(f"Output: {output_file}")
    
    try:
        # Process the point cloud
        las = fPCpretreatment(
            chunk=str(input_path),
            classify=args.classify,
            LMA=LMA,
            WD=WD,
            WD_bush=args.WD_bush,
            LMA_bush=args.LMA_bush,
            H_strata_bush=args.H_strata_bush,
            Height_filter=args.Height_filter,
            start_date=args.start_date,
            season_filter=args.season_filter,
            deviation_days=args.deviation_days,
            plot_hist_days=args.plot_hist_days
        )
        
        if las is None:
            print("Warning: Processing returned None. No output file created.", file=sys.stderr)
            sys.exit(1)
        
        # Write output file
        las.write(str(output_path))
        print(f"Successfully processed and saved to: {output_file}")
        
    except Exception as e:
        print(f"Error during processing: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

