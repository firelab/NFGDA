import numpy as np
from osgeo import gdal, osr


# parser = argparse.ArgumentParser(
#     description="Convert radar NPZ data to GeoTIFF using radar lat/lon."
# )

# parser.add_argument(
#     "npz_path",
#     type=str,
#     help="Path to input NPZ file"
# )

# parser.add_argument(
#     "radar_lat",
#     type=float,
#     help="Radar latitude in degrees"
# )

# parser.add_argument(
#     "radar_lon",
#     type=float,
#     help="Radar longitude in degrees"
# )

# args = parser.parse_args()

# npz_path = args.npz_path
# radar_lat = args.radar_lat
# radar_lon = args.radar_lon

# #------------example input data for testing:
# #npz_path = "nf_predKABX20200707_012805_V06.npz"
# #KABX radar location
# #radar_lat = 35.149722    # degrees
# #radar_lon = -106.823889  # degrees

# # ---------------------------
# # Other parameters
# # ---------------------------
# utm_tif = "radar_reflectivity_utm.tif"
# final_tif = "radar_reflectivity_latlon.tif"

# pixel_size_m = 500.0   # 500 m spacing
# channel_index = 1      # channel 1 (0-based)

# # ---------------------------
# # Load data
# # ---------------------------
# data = np.load(npz_path)
# array = data['inputNF'] 

# # Flip vertically (upside down)
# array = np.flipud(array)

# refl = array[:, :, channel_index].astype(np.float64)
# ny, nx = refl.shape
# print(f"ny, nx = {ny, nx}")

# ---------------------------
# Determine best-fit UTM zone
# ---------------------------

class Gdal_Writer:
    def __init__(self, radar_lat, radar_lon, pixel_size_m, nx, ny):
        utm_zone = int((radar_lon + 180) / 6) + 1
        is_northern = radar_lat >= 0
        epsg_utm = 32600 + utm_zone if is_northern else 32700 + utm_zone

        utm_srs = osr.SpatialReference()
        utm_srs.ImportFromEPSG(epsg_utm)

        wgs84 = osr.SpatialReference()
        wgs84.ImportFromEPSG(4326)

        to_utm = osr.CoordinateTransformation(wgs84, utm_srs)
        center_x, center_y, _ = to_utm.TransformPoint(radar_lat, radar_lon)
        origin_x = center_x - (nx / 2) * pixel_size_m
        origin_y = center_y + (ny / 2) * pixel_size_m
        self.geotransform = (
            origin_x,
            pixel_size_m,
            0.0,
            origin_y,
            0.0,
            -pixel_size_m
        )
        self.utm_wkt = utm_srs.ExportToWkt()

    def log_geo_tif(self, fn, array):
        mem_driver = gdal.GetDriverByName("MEM")

        utm_ds = mem_driver.Create(
            "", # No filename needed for MEM driver
            array.shape[1], # nfgda var [y,x,v] in here order for [nx,ny,nz]=[1,0,2]
            array.shape[0], 
            array.shape[2], 
            gdal.GDT_Float64
        )
        utm_ds.SetGeoTransform(self.geotransform)
        utm_ds.SetProjection(self.utm_wkt)
        for iv in range(array.shape[-1]):
            band = utm_ds.GetRasterBand(iv+1) ## band index start from 1
            band.WriteArray(np.flipud(array[:, :, iv]).astype(np.float64))
            band.SetNoDataValue(-9999)
        gdal.Warp(
            fn,
            utm_ds,
            dstSRS="EPSG:4326",
            srcNodata=-9999, dstNodata=-9999,
            resampleAlg=gdal.GRA_NearestNeighbour,
            format="GTiff",
            creationOptions=["COMPRESS=LZW", "TILED=YES"]
        )
        utm_ds = None

# wt = Gdal_Writer( radar_lat, radar_lon, pixel_size_m, nx, ny)
# wt.log_geo_tif(final_tif,array)
