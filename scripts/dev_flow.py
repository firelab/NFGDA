import asyncio
from concurrent.futures import ProcessPoolExecutor
import sys
import datetime
import os
import numpy as np
import nfgda.NF_Lib as NF_Lib
from nfgda import nf_path
from nfgda.NFGDA_load_config import *

import nfgda.NF_Lib as NF_Lib
import matplotlib.pyplot as plt

last_nexrad = custom_start_time-datetime.timedelta(minutes=1)
# exit_time = custom_end_time

scans = NF_Lib.aws_int.get_avail_scans_in_range(last_nexrad, last_nexrad+datetime.timedelta(minutes=20), radar_id)
# for vol in scans:
#     NF_Lib.get_nexrad (path_config,vol)
# for iv in range(len(scans)-1):
#     NF_Lib.nfgda_unit_step(scans[iv].filename,scans[iv+1].filename)

# for iv in range(len(scans)-2):
#     NF_Lib.nfgda_forecast(scans[iv+1].filename,scans[iv+2].filename)


# filepath = nf_path.get_nf_forecast_npz_name(scans[1].filename, path_config)

# conn = NF_Lib.Prediction_Connection.load(filepath)

# st = conn.timestamp.astype('datetime64[m]')
# second_offset = (st.astype(int)*60) % forecast_step_sec
# tvec = st - second_offset*np.timedelta64(1, 's') \
#     +np.arange(forecast_step_sec,forecast_period_sec+1,forecast_step_sec)*np.timedelta64(1, 's')

# print(st)
# print(tvec)
conns = []
for iv in range(len(scans)-2):
    filepath = nf_path.get_nf_forecast_npz_name(scans[iv+1].filename, path_config)
    conns.append(NF_Lib.Prediction_Connection.load(filepath))
    print(conns[-1].igp_anchor.shape)

arcs = NF_Lib.get_forecast(conns[-1],conns[-1].timestamp+np.timedelta64(180, 's'))
arcs[0].write_geojson('test.geojson')
# NF_Lib.nfgda_stochastic_summary(conns,scans[iv+1].filename)