import os
import datetime

def get_nf_input_name(l2_file, path_config):
    fn = l2_file.split('.')[0]+'.npz'
    return os.path.join(path_config.nf_dir, fn)

def get_nf_detection_name(l2_file, path_config,ext='npz'):
    fn = 'nf_pred'+l2_file+'.'+ext
    return os.path.join(path_config.nf_preds_dir, fn)

def get_nf_forecast_name(l2_file, path_config, valid_time):
    savedir = os.path.join(path_config.nf_forecast_dir, 'forecast-'+l2_file)
    os.makedirs(savedir,exist_ok=True)
    fn = f'NFGDA-forecast-{l2_file[:4]}'+valid_time.astype(datetime.datetime).strftime('%Y%m%d_%H%M%S')+'.png'
    return os.path.join(savedir, fn)

def get_nf_forecast_pkl_name(l2_file, path_config):
    savedir = os.path.join(path_config.nf_forecast_dir, 'forecast-'+l2_file)
    os.makedirs(savedir,exist_ok=True)
    fn = f'NFGDA-forecast-{l2_file}.pkl'
    return os.path.join(savedir, fn)

def get_nf_forecast_npz_name(l2_file, path_config):
    savedir = os.path.join(path_config.nf_forecast_dir)
    # os.makedirs(savedir,exist_ok=True)
    fn = f'NFGDA-forecast-{l2_file}.npz'
    return os.path.join(savedir, fn)

def get_nf_s_forecast_name(path_config, valid_time, ext='png'):
    savedir = os.path.join(path_config.nf_forecast_dir, 'forecast-summary')
    fn = f'NFGDA-forecast-{path_config.radar_id}'+valid_time.astype(datetime.datetime).strftime('%Y%m%d_%H%M%S')+'.'+ext
    return os.path.join(savedir, fn)