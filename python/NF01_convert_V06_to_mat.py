import os
import pyart
import numpy as np
import argparse
import scipy
import configparser

# def read_info_from_radar_name(radar_file):
#     """ Copied from https://github.com/PreciousJatau47/VAD_correction/blob/master/RadarHCAUtils.py
#     Finds the corresponding hca PPI name given a radar PPI name.
#     :param radar_file:
#     :return:
#     """
#     radar_name = radar_file[:4]
#     year = radar_file[4:8]
#     month = radar_file[8:10]
#     day = radar_file[10:12]
#     hh = radar_file[13:15]
#     mm = radar_file[15:17]
#     ss = radar_file[17:19]
#     return radar_name, year, month, day, hh, mm, ss

def ReadSliceElevation(radar, slice_idx):
    """ Copied from https://github.com/PreciousJatau47/VAD_correction/blob/master/RadarHCAUtils.py
    :param radar:
    :param slice_idx:
    :return:
    """
    sweep_ind = radar.get_slice(slice_idx)
    radar_el = radar.elevation['data'][sweep_ind]
    return radar_el


def ReadRadarSliceUpdate(radar, slice_idx):
    """ Copied from https://github.com/PreciousJatau47/VAD_correction/blob/master/RadarHCAUtils.py
    :param radar:
    :param slice_idx:
    :return:
    """
    radar_range = radar.range['data'] / 1000  # in km
    sweep_ind = radar.get_slice(slice_idx)
    radar_az_deg = radar.azimuth['data'][sweep_ind]  # in degrees
    radar_el = radar.elevation['data'][sweep_ind]

    ref_shape = radar.fields["reflectivity"]['data'][sweep_ind].shape
    placeholder_matrix = np.full(ref_shape, np.nan, dtype=np.float64)
    placeholder_mask = np.full(ref_shape, False, dtype=bool)

    data_slice = []
    labels_slice = list(radar.fields.keys())
    labels_slice.sort()
    mask_slice = []
    var_mask_slice = []

    for radar_product in labels_slice:
        if np.sum(radar.fields[radar_product]['data'][sweep_ind].mask == False) > 0:
            data_slice.append(radar.fields[radar_product]['data'][sweep_ind])
            mask_slice.append(True)
            var_mask_slice.append(radar.fields[radar_product]['data'][sweep_ind].mask)
        else:
            data_slice.append(placeholder_matrix)
            mask_slice.append(False)
            var_mask_slice.append(placeholder_mask)

    return radar_range, radar_az_deg, radar_el, data_slice.copy(), mask_slice.copy(), labels_slice, var_mask_slice


END_GATE = 400
NUM_AZ = 720


def convert_v06_to_mat(v06_folder, case_id, mat_folder, i_start, i_end):
    var_2_parrot_idx = {'reflectivity': 0, 'velocity': 1, 'spectrum_width': 2, 'differential_phase': 3,
                        'cross_correlation_ratio': 4, 'differential_reflectivity': 5}
    # target_el = 1.25
    v06_folder = os.path.join(v06_folder, case_id)
    l2_files = [entry for entry in os.listdir(v06_folder) if entry.endswith("V06")]
    l2_files.sort()

    for i in range(len(l2_files)):
        l2_file = l2_files[i]
        if not (l2_file.endswith('_V06') or l2_file.startswith('._')):
            continue

        if i < i_start or i > i_end:
            continue

        print("case number: ", i)
        print("Processing ", l2_file)

        # Output path.
        output_folder = os.path.join(mat_folder, 'POLAR', case_id)

        if not os.path.isdir(output_folder):
            os.makedirs(output_folder)

        # mat_file = "".join(['polar', case_id, '%02i' % (i + 1), '.mat'])
        mat_file = f'polar_{i + 1:0>2d}_' + l2_file.split('.')[0]+'.mat'

        output_path = os.path.join(output_folder, mat_file)

        py_path = os.path.join(output_folder, mat_file[:-4]+'.npz')

        # read l2 data
        radar_obj = pyart.io.read_nexrad_archive(os.path.join(v06_folder, l2_file))

        # TODO(pjatau) erase below.
        nsweeps = radar_obj.nsweeps
        vcp = radar_obj.metadata['vcp_pattern']
        print("VCP: ", vcp)

        # VCP 212.
        # slices 0-2-4 contain only dual-pol. super res.
        # slices 1-3-5 contain vel products. super res.
        # slices >= 6 contain all products. normal res.

        # Initialize data cube
        # PARROT = np.ma.full((END_GATE, NUM_AZ, 6), np.nan, dtype=np.float64)
        PARROT = np.ma.array(np.ma.array(np.full((END_GATE, NUM_AZ, 6), np.nan, dtype=np.float64)), mask=np.full((END_GATE, NUM_AZ, 6), True))
        in_parrot = np.full(6,False)

        for slice_idx in range(nsweeps):
            radar_el = ReadSliceElevation( radar_obj, slice_idx)
            scan_el = np.nanmedian(radar_el)
            # if abs(scan_el-target_el)>0.3:
            #     continue
            radar_range, az_sweep_deg, radar_el, data_slice, mask_slice, labels_slice, data_mask_slice = ReadRadarSliceUpdate(
                radar_obj, slice_idx)
            print("Processing elevation {} degrees".format(np.nanmedian(radar_el)))

            i_zero_az = np.argmin(np.abs(az_sweep_deg))
            az_shift = -i_zero_az

            var_idx_slice = {labels_slice[i]: i for i in range(len(labels_slice))}

            for var in var_2_parrot_idx.keys():
                i_var = var_idx_slice[var]
                i_parrot = var_2_parrot_idx[var]

                if not mask_slice[i_var] or in_parrot[i_parrot]:
                    continue
                in_parrot[i_parrot] = True

                print("Processing {}. parrot idx {}".format(var, i_parrot))

                curr_data = data_slice[i_var][:, :END_GATE]
                curr_mask = data_mask_slice[i_var][:, :END_GATE]
                # curr_data[curr_mask] = np.nan  # (720, 400)
                curr_data = np.roll(a=curr_data, shift=az_shift, axis=0)
                PARROT[:, :, i_parrot] = curr_data.T
                # PARROT[:, :, i_parrot].data = curr_data.T
                # PARROT[:, :, i_parrot].mask = curr_mask.T
                # print(PARROT.mask)
                # exit()
            if np.min(in_parrot):
                timestamp=np.datetime64(pyart.graph.common.generate_radar_time_sweep(radar_obj,slice_idx))
                print("slice idx {} timestamp".format(slice_idx),timestamp)
                break
            print()
        print()
        # scipy.io.savemat(output_path, {"PARROT": PARROT})
        np.savez(py_path, PARROT=PARROT.data,mask=PARROT.mask,timestamp=timestamp)

    return


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description="")
    # parser.add_argument('-v06_folder', type=str, default="../V06", help='folder containing V06 files')
    # parser.add_argument('-case_id', type=str, help='V06 case folder')
    # parser.add_argument('-mat_folder', type=str, default="../mat", help='output folder for mat files')
    # parser.add_argument('-i_start', type=int, default=0, help='0-based index of first V06 file')
    # parser.add_argument('-i_end', type=int, default=99, help='0-based index of last V06 file')
    # args = parser.parse_args()
    # convert_v06_to_mat(v06_folder=args.v06_folder, case_id=args.case_id, mat_folder=args.mat_folder,
    #                    i_start=args.i_start, i_end=args.i_end)
    config = configparser.ConfigParser()
    config.read("./NFGDA.ini")
    case_name = config["Settings"]["case_name"]
    convert_v06_to_mat(v06_folder="../V06", case_id=case_name, mat_folder="../mat",
                       i_start=0, i_end=99)
