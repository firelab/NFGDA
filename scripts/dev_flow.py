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
for vol in scans:
    NF_Lib.get_nexrad (path_config,vol)
for iv in range(len(scans)-1):
    NF_Lib.nfgda_unit_step(scans[iv].filename,scans[iv+1].filename)

for iv in range(len(scans)-2):
    NF_Lib.nfgda_forecast(scans[iv+1].filename,scans[iv+2].filename)

# worker.save_conn(0,nf_path.get_nf_forecast_npz_name(l2_file_0, path_config))
# l2_file_0
# path_config
filepath = nf_path.get_nf_forecast_npz_name(scans[1].filename, path_config)

conn = NF_Lib.Prediction_Connection.load(filepath)

st = conn.timestamp.astype('datetime64[m]')
second_offset = (st.astype(int)*60) % forecast_step_sec
tvec = st - second_offset*np.timedelta64(1, 's') \
    +np.arange(forecast_step_sec,forecast_period_sec+1,forecast_step_sec)*np.timedelta64(1, 's')

print(st)
print(tvec)
conns = []
for iv in range(len(scans)-2):
    filepath = nf_path.get_nf_forecast_npz_name(scans[iv+1].filename, path_config)
    conns.append(NF_Lib.Prediction_Connection.load(filepath))

# py_path = nf_path.get_nf_detection_name(scans[iv+1].filename, path_config)
# data = np.load(py_path)
# tnow = st +( -second_offset + forecast_step_sec )*np.timedelta64(1, 's')
# pdata = np.ma.masked_where(rmask,data['inputNF'][:,:,1])

# for t in tvec:
#     fig, axs = plt.subplots(1, 1, figsize=(3.3/0.7, 3/0.7),dpi=150)
#     pcz=axs.pcolormesh(Cx,Cy,pdata,cmap=cl.zmap,norm=cl.znorm)
#     axs.set_xlim(-100,100)
#     axs.set_ylim(-100,100)
#     axs.set_xlabel('x(km)')
#     axs.set_ylabel('y(km)',labelpad=-10)
#     axs.set_aspect('equal')
#     forecast_anchors = []
#     for conn in conns:
#         if t < conn.timestamp:
#             print(f'{t} < {conn.timestamp}')
#             continue
#         anchors = NF_Lib.get_forecast(conn,t)
#         axs.plot(anchors[0].arc_anchors[:,0,:].T,anchors[0].arc_anchors[:,1,:].T,alpha=0.7,color='k')
#         axs.plot(anchors[1].arc_anchors[:,0,:].T,anchors[1].arc_anchors[:,1,:].T,alpha=0.7,color='r')
#         # fig.suptitle(t.astype(datetime.datetime).strftime('%Y/%m/%d %H:%M:%S')+f' (+{int(dt)} mins)',y=0.97)
#         dt = (t-conn.timestamp)/np.timedelta64(60, 's')
#         fig.suptitle(tnow.astype(datetime.datetime).strftime('%Y/%m/%d %H:%M:%S')+'\n'
#             +t.astype(datetime.datetime).strftime('%Y/%m/%d %H:%M:%S')+f' (+{int(dt)} mins)',y=0.97)
#     # savedir = os.path.join(path_config.nf_forecast_dir)
#     # os.makedirs(savedir,exist_ok=True)
#         fn = f'NFGDA-forecast-'+t.astype(datetime.datetime).strftime('%Y%m%d_%H%M%S')+'.png'
#         fig.savefig(os.path.join(path_config.nf_forecast_dir, fn))
#         for ln in axs.lines[:]:
#             ln.remove()
#     plt.close(fig)

NF_Lib.nfgda_stochastic_summary(conns,scans[iv+1].filename)