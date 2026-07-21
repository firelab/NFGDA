import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
import nfgda.NF_Lib as NF_Lib
from pyproj import Transformer
from nfgda import colorlevel as cl
from matplotlib.colors import ListedColormap
import os
from datetime import datetime

# 1. Open the raster with GDAL
Detection_path = "../plot_res/nf_predKABX20200707_013332_V06.tif"
forecast_path = "../plot_res/NFGDA-forecast-KABX20200707_014700.tif"
ds_det = gdal.Open(Detection_path)
gt = ds_det.GetGeoTransform()

width = ds_det.RasterXSize
height = ds_det.RasterYSize

# Build cell coordinate arrays (Cx, Cy)
x_coords = gt[0] + np.arange(width + 1) * gt[1]
y_coords = gt[3] + np.arange(height + 1) * gt[5]
Cx, Cy = np.meshgrid(x_coords, y_coords)

# Read Band 2 for reflectivity shading
band1 = ds_det.GetRasterBand(2)
bg_data = band1.ReadAsArray()
nodata = band1.GetNoDataValue()
bg_data = np.ma.masked_equal(bg_data, nodata)

# Read Band 7 for NFGDA detection overlay (val == 1 will be red)
band7 = ds_det.GetRasterBand(7)
bin_data = band7.ReadAsArray()
bin_data[bin_data < 0] = 0

# -------------------------------------------------------------
# 2. Load Forecast Proxy Raster (Band 1 for blue contours)
# -------------------------------------------------------------
ds_fcst = gdal.Open(forecast_path)
fcst_data = ds_fcst.GetRasterBand(1).ReadAsArray()

gt_fcst = ds_fcst.GetGeoTransform()
fcst_x = gt_fcst[0] + np.arange(ds_fcst.RasterXSize + 1) * gt_fcst[1]
fcst_y = gt_fcst[3] + np.arange(ds_fcst.RasterYSize + 1) * gt_fcst[5]
FCx, FCy = np.meshgrid(fcst_x, fcst_y)

# -------------------------------------------------------------
# 3. Set up the Plot colormap/norm
# -------------------------------------------------------------
fig, axs = plt.subplots(1, 1, figsize=(3.3/0.7, 3/0.7), dpi=150)

# Background shading
pcz = axs.pcolormesh(Cx, Cy, bg_data, cmap=cl.zmap, norm=cl.znorm, shading='flat')

# Overlay: Detection Band 7 binary pixels where val == 1 displayed as red
bin_masked = np.ma.masked_where(bin_data != 1, bin_data)
red_cmap = ListedColormap(['red'])
axs.pcolormesh(Cx, Cy, bin_masked, cmap=red_cmap, shading='flat', alpha=0.9)

# Overlay: Forecast Band 1 as blue contours on top
axs.contour(FCx[:-1, :-1], FCy[:-1, :-1], fcst_data, levels=[30], colors='blue', linewidths=1.5)

axs.set_xlabel("Lon")
axs.set_ylabel("Lat")

# parse title
station = Detection_path.split('_')[-3][4:4+4]
det_date = Detection_path.split('_')[-3][8:8+8]
det_time = Detection_path.split('_')[-2]
det_dt = datetime.strptime(f"{det_date} {det_time}", "%Y%m%d %H%M%S")

fcst_date = forecast_path.split('-')[-1][4:4+8]
fcst_time = forecast_path.split('-')[-1].split('_')[-1][:6]
fcst_dt = datetime.strptime(f"{fcst_date} {fcst_time}", "%Y%m%d %H%M%S")

# Format for title
det_str = det_dt.strftime("%Y-%m-%d %H:%M:%S")
fcst_str = fcst_dt.strftime("%Y-%m-%d %H:%M:%S")

plot_title = f"{station} | Detection : {det_str} \n-> Forecast : {fcst_str} UTC"
axs.set_title(plot_title)

plt.show()
