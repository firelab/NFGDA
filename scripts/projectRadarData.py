#!/usr/bin/env python3

import argparse
import numpy as np
import nfgda.nf_gdal

parser = argparse.ArgumentParser(
    description="Convert radar NPZ data to GeoTIFF using radar lat/lon."
)

parser.add_argument(
    "npz_path",
    type=str,
    help="Path to input NPZ file"
)

parser.add_argument(
    "radar_lat",
    type=float,
    help="Radar latitude in degrees"
)

parser.add_argument(
    "radar_lon",
    type=float,
    help="Radar longitude in degrees"
)

args = parser.parse_args()

npz_path = args.npz_path
radar_lat = args.radar_lat
radar_lon = args.radar_lon
"""
Example CLI for testing:
    $ ./projectRadarData.py nf_predKABX20200707_012805_V06.npz 35.149722 -106.823889
"""
final_tif = "radar_latlon.tif"

pixel_size_m = 500.0   # 500 m spacing

data = np.load(npz_path)
array = data['inputNF']

wt = nfgda.nf_gdal.Gdal_Writer( radar_lat, radar_lon, pixel_size_m, array.shape[1], array.shape[0])
wt.log_geo_tif(final_tif,array)