import configparser
import numpy as np
import os
import datetime
from scipy.interpolate import LinearNDInterpolator
from importlib.resources import files
from . import math_kit as mk
from . import colorlevel as cl
VM=cl.VarMap()
varname_table=VM.varname_table
varunit_table=VM.varunit_table

config = configparser.ConfigParser()
config.read("NFGDA.ini")
radar_id = config["Settings"]["radar_id"]
radar_lat, radar_lon = list(map(float,config.get("labels", "loc").split(",")))
export_preds_dir = config["Settings"]["export_preds_dir"]
evalbox_on = config.getboolean('Settings', 'evalbox_on')
export_forecast_dir = config["Settings"]["export_forecast_dir"]
V06_dir = config["Settings"]["V06_dir"]
custom_start_time = config["Settings"]["custom_start_time"]
custom_end_time = config["Settings"]["custom_end_time"]
use_gdal = True
if use_gdal:
    from . import nf_gdal
if len(custom_start_time.split(','))==6:
    custom_start_time = datetime.datetime(*map(int,custom_start_time.split(',')), tzinfo=datetime.timezone.utc)
    custom_end_time = datetime.datetime(*map(int,custom_end_time.split(',')), tzinfo=datetime.timezone.utc)
else:
    custom_start_time = None
    custom_end_time = None

# fig_dir = config["Settings"]["fig_dir"]
# label_on = config.getboolean('labels', 'label_on')
# if label_on:
    
#     radar_loc = list(map(float,config.get("labels", "rloc").split(",")))
#     sitex, sitey = mk.geopoints_to_relative_xy(radar_loc,label_loc)
pixel_size_m = 500
Cx, Cy = np.meshgrid(np.arange(-100,100.5,pixel_size_m/1e3),np.arange(-100,100.5,pixel_size_m/1e3))
r = np.sqrt(Cx**2+Cy**2)
rmask = r>=100

os.makedirs(V06_dir,exist_ok=True)
nf_dir = V06_dir+'npz/'
os.makedirs(nf_dir,exist_ok=True)
os.makedirs(export_preds_dir,exist_ok=True)
sf_dir = os.path.join(export_forecast_dir, 'forecast-summary')
os.makedirs(sf_dir,exist_ok=True)

PARROT_mask_on = False

thrREF = -5
thrdREF = 0.3
RegR = np.arange(0,400)/4
RegAZ = np.arange(0,360,0.5)*np.pi/180
RegPolarX = RegR[:,np.newaxis] * np.sin(RegAZ[np.newaxis,:])
RegPolarY = RegR[:,np.newaxis] * np.cos(RegAZ[np.newaxis,:])

###### Beta Cell magic numbers ######
cellthresh = 5
cbcellthrsh = 0.8
cellcsrthresh=0.5
crsize = 5
cellINT = crsize + 2
widecellINT =crsize+4
avgINT = 8
s2xnum = [10, 15]
s2ynum = [-3, 1]
s2xdel = s2xnum[1]-s2xnum[0]
s2ydel = s2ynum[1]-s2ynum[0]
s2g = s2ydel/s2xdel
s2gc = s2ynum[1]-s2g*s2xnum[1]
Celldp = np.load(files("nfgda").joinpath("Celldp.npy"))
Celldpw = np.load(files("nfgda").joinpath("Celldpw.npy"))
###### Beta Cell magic numbers ######

###### FTC Beta Z, dZ displacements ######
datacy = np.arange(-8,9).reshape(1,-1)
datacx = np.zeros((1,17))
datac = np.swapaxes(np.array([datacy,datacx]),0,2)

datasy = np.array([*np.arange(-7,0,2),0,*np.arange(1,8,2),*np.arange(-7,0,2),0,*np.arange(1,8,2)]).reshape(1,-1)
datasx = np.array([-4*np.ones((9)),4*np.ones((9))]).reshape(1,-1)
datas = np.swapaxes(np.array([datasy,datasx]),0,2)
###### FTC Beta Z, dZ displacements ######

class path_struct():
    def __init__(self):
        self.radar_id = radar_id
        self.nf_dir = nf_dir
        self.V06_dir = V06_dir
        self.nf_preds_dir = export_preds_dir
        self.nf_forecast_dir = export_forecast_dir
        if use_gdal:
            self.gdal_writer = nf_gdal.Gdal_Writer( radar_lat, radar_lon, pixel_size_m, Cx.shape[1], Cx.shape[0])

path_config = path_struct()

###### Utilities ######
def tprint(*args, **kwargs):
    print(f"[{datetime.datetime.now():%H:%M:%S}]", *args, **kwargs)

class C:
    RESET   = "\033[0m"
    # Standard colors
    BLACK   = "\033[30m"
    RED     = "\033[31m"
    GREEN   = "\033[32m"
    YELLOW  = "\033[33m"
    BLUE    = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN    = "\033[36m"
    WHITE   = "\033[37m"

    # Bright (bold) colors
    BLACK_B   = "\033[1;30m"
    RED_B     = "\033[1;31m"
    GREEN_B   = "\033[1;32m"
    YELLOW_B  = "\033[1;33m"
    BLUE_B    = "\033[1;34m"
    MAGENTA_B = "\033[1;35m"
    CYAN_B    = "\033[1;36m"
    WHITE_B   = "\033[1;37m"
###### Utilities ######
dl_tag = f"{C.CYAN}[Downloader]{C.RESET} "
cv_tag = f"{C.BLUE}[Converter]{C.RESET} "
ng_tag = f"{C.GREEN}[NFGDA]{C.RESET} "
df_tag = f"{C.YELLOW}[FORECAST]{C.RESET} "
sf_tag = f"{C.YELLOW_B}[S FORECAST]{C.RESET} "
ht_tag = f"{C.MAGENTA}[Host]{C.RESET} "
