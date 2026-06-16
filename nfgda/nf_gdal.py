import numpy as np
from osgeo import gdal, osr

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
        if array.ndim ==3:
            Nvar = array.shape[2]
        else:
            Nvar = 1
        utm_ds = mem_driver.Create(
            "", # No filename needed for MEM driver
            array.shape[1], # nfgda var [x,y,var] convert to GeoTIFF Raster Image [vertical(NS),horizontal(EW),band] // [0,1,2] --> [1,0,2]
            array.shape[0], 
            Nvar, 
            gdal.GDT_Float64
        )
        utm_ds.SetGeoTransform(self.geotransform)
        utm_ds.SetProjection(self.utm_wkt)
        if Nvar != 1:
            for iv in range(array.shape[-1]):
                band = utm_ds.GetRasterBand(iv+1) ## band index start from 1
                band.WriteArray(np.flipud(array[:, :, iv]).astype(np.float64))
                band.SetNoDataValue(-9999)
        else:
            band = utm_ds.GetRasterBand(1) ## band index start from 1
            band.WriteArray(np.flipud(array).astype(np.float64))
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
