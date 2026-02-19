from numpy.lib.stride_tricks import sliding_window_view
import scipy.io
from scipy.signal import medfilt2d
from scipy.spatial.distance import pdist, squareform
from scipy.ndimage import convolve
# from scipy.ndimage import gaussian_filter
from skimage.morphology import skeletonize, disk, binary_dilation, remove_small_objects
from skimage.measure import label
import matplotlib.pyplot as plt
import sys
import pyart
from pathlib import Path
from . import nf_path
from .NFGDA_load_config import *
import pickle
if use_gdal:
    from . import nf_gdal

import nexradaws
aws_int = nexradaws.NexradAwsInterface()
interpolator = LinearNDInterpolator((RegPolarX.reshape(-1),RegPolarY.reshape(-1)), np.zeros(RegPolarX.shape).reshape(-1))

kernel = np.ones((3,3), dtype=int)
ele_t_const, mean_t_const = 29.04, 28.43

def rot_displace(dp,origindeg):
    dpvector = np.swapaxes(dp,1,2)
    origindeg = origindeg*np.pi/180
    backprocess = np.array([[np.cos(origindeg),np.sin(origindeg)], \
                            [-np.sin(origindeg),np.cos(origindeg)]])
    rotcord = np.matmul(backprocess,dpvector)
    rotidx = np.round(rotcord)
    return np.swapaxes(rotidx,1,2)

ftcc=[]
ftcs=[]

rotdegree = 180/9
for irot in np.arange(0,180,rotdegree):
    ftcc.append(rot_displace(datac,irot))
    ftcs.append(rot_displace(datas,irot))
ftcc = np.array(ftcc)
ftcs = np.array(ftcs)
######### FTC Beta Z, dZ displacements ###########
mvdiscx = np.zeros((17,17))
mvdiscy = np.zeros((17,17))
for ix in range(17):
    mvdiscx[ix,:]=np.ceil(np.arange(-8,9)*np.sin(np.pi/2/8*(ix)))
    mvdiscy[ix,:]=np.ceil(np.arange(-8,9)*np.cos(np.pi/2/8*(ix)))
nccx = mvdiscx.shape[0]
mvdisc = np.swapaxes(np.array([mvdiscy,mvdiscx]),0,2)[:,:,np.newaxis,:]

def make_ftc_cscore(c_para):
    cnum1, cnum2, csig1, cfactor1, cintersec1, csig2, cfactor2, cintersec2, cyfill = c_para
    def f(cbox):
        # params is captured from outer scope
        llscore = np.zeros(cbox.shape)
        # llscore = np.full(cbox.shape,np.nan)
        llscore[cbox<=cnum1] = gaussmf(cbox[cbox<=cnum1], csig1, cnum1)*cfactor1+cintersec1
        llscore[np.logical_and(cbox>cnum1, cbox<=cnum2)] = cyfill
        llscore[cbox>cnum2] = gaussmf(cbox[cbox>cnum2], csig2, cnum2)*cfactor2+cintersec2
        return llscore
    return f

def make_ftc_sscore(s_para):
    snum1, snum2, ssig1, sfactor1, sintersec1, ssig2, sfactor2, sintersec2, syfill = s_para
    def f(sbox):
        # params is captured from outer scope
        ssscore = np.zeros(sbox.shape)
        # ssscore = np.full(sbox.shape,np.nan)
        ssscore[sbox<snum1] = syfill
        con1 = np.logical_and(sbox>=snum1, sbox<=snum2)
        con2 = sbox>snum2
        ssscore[con1] = gaussmf(sbox[con1], ssig1, snum1)*sfactor1 + sintersec1
        ssscore[con2] = gaussmf(sbox[con2], ssig2, snum2)*sfactor2 + sintersec2
        return ssscore
    return f

class FTC_PLAN:
    def __init__(self,displace,scorefun,scale):
        self.displace = displace
        self.scorefun = scorefun
        self.numINT = np.max(np.abs(displace))
        self.scale = scale
    def gather_pixels(self,ar,center):
        # idx [direction, displacement, center_pixel, yx]
        idx = (center[np.newaxis,np.newaxis,:,:] + self.displace).astype(int)
        return ar[idx[:,:,:,0],idx[:,:,:,1]]
    def get_score(self,ar,center):
        cbox = self.gather_pixels(ar,center)
        pixel_score = self.scorefun(cbox)
        # pixel_score[np.isnan(pixel_score)] = -3
        return np.nansum(pixel_score,axis=1)

def gen_beta(a2,a2_thr,ftcs):
    center_indices = np.argwhere(a2>a2_thr)
    c_indices = clean_indices(center_indices, a2.shape, ftcs[0].numINT)
    total_score = ftcs[0].get_score(a2,c_indices)
    score_scale = ftcs[0].scale*ftcs[0].displace.shape[1]
    for ftcplan in ftcs[1:]:
        total_score += ftcplan.get_score(a2,c_indices)
        score_scale += ftcplan.scale*ftcplan.displace.shape[1]
    total_score = np.max(total_score,axis=0)
    scoremt = np.zeros(a2.shape)
    scoremt[c_indices[:,0],c_indices[:,1]] = total_score/score_scale
    return scoremt

def gaussmf(x, sigma, c):
    return np.exp(-((x - c) ** 2) / (2 * sigma ** 2))

def probor(ar):
    buf = np.zeros(ar.shape[:-1])
    for iv in range(ar.shape[-1]):
        buf = buf + ar[...,iv] - buf * ar[...,iv]
    return buf

def clean_indices(idx,shp,edg):
    dim0 = idx[:,0]
    dim1 = idx[:,1]
    inbox = (dim0>=edg) & (dim0< shp[0]-edg) & (dim1>=edg) & (dim1< shp[1]-edg)
    return idx[inbox,:]

def post_moving_avg(a2):
    center_indices = np.argwhere(a2>0)

    c_indices = clean_indices(center_indices, a2.shape, avgINT)
    cidx = (c_indices[np.newaxis,np.newaxis,...] + mvdisc).astype(int)
    cbox = a2[cidx[...,0],cidx[...,1]]
    # cbr = np.sum(cbox>thrREF,0)/datacx.size
    # cbox = cbox[:,cbr>0.5]
    cbr = np.sum(cbox>0,axis=0)/nccx
    validcenter = np.max(cbr>0.1,axis=0)
    cbox = cbox[:,:,validcenter,...]
    mc = np.nanmean(cbox,axis=0)
    mc[np.logical_not(cbr>0.1)]=0
    s_indices = c_indices[validcenter,:]
    result = np.zeros(a2.shape)
    result[s_indices[:,0],s_indices[:,1]] = np.max(mc,axis=0)
    return result

class NFModule:
    def __init__(self,fismat):
        buf = scipy.io.loadmat(fismat)
        # [1,1,rule,vars]
        self.c = buf['incoef'][:,1,:][np.newaxis,np.newaxis,...]
        self.sig = buf['incoef'][:,0,:][np.newaxis,np.newaxis,...]
        self.outcoef = buf['outcoef'][np.newaxis,np.newaxis,...]
        self.rulew = buf['rulelogic'][:,0]
        self.rulecon = buf['rulelogic'][:,1]
    def eval_fis(self,pxls):
        # pxls[x,y,vars] -> [x,y,rule,vars]
        x = pxls[:,:,np.newaxis,:]
        irr = np.exp(-(x-self.c)**2/(2*self.sig**2))
        # w[x,y,rule]
        w = np.zeros(irr.shape[:-1])
        for ir in range(self.rulecon.size):
            if self.rulecon[ir]==1:
                w[...,ir] = np.prod(irr[...,ir,:],axis=-1)
            else:
                w[...,ir] = probor(irr[...,ir,:])
        sw = np.sum(w,axis=-1)
        pad_window = tuple([(0,0),] * (x.ndim-1) + [(0,1),])
        orr = np.pad(x, pad_width = pad_window, mode='constant', constant_values=1)
        # orr = [x,ones(size(x,1),1)];
        unImpSugRuleOut = np.sum(orr*self.outcoef,axis=-1)
        return np.sum( unImpSugRuleOut*w, axis=-1 )/sw

fuzzGST = NFModule(files("nfgda").joinpath("NF00ref_YHWANG_fis4python.mat"))

def nfgda_unit_step(l2_file_0,l2_file_1,process_tag=''):
    ifn = l2_file_1
    tprint(ng_tag+f'{l2_file_0} -> {l2_file_1}')

    buf = np.load(nf_path.get_nf_input_name(l2_file_0,path_config))
    PARROT0 = buf['PARROT']
    if PARROT_mask_on:
        PARROT0[buf['mask']] = np.nan
    PARROT0 = np.asfortranarray(PARROT0)

    PARROT_buf = np.load(nf_path.get_nf_input_name(l2_file_1,path_config))
    PARROT = PARROT_buf['PARROT']
    if PARROT_mask_on:
        PARROT[PARROT_buf['mask']] = np.nan
    PARROT = np.asfortranarray(PARROT)
    
    diffz = PARROT[:,:,0] - PARROT0[:,:,0]
    
    PARITP = np.zeros((*Cx.shape,PARROT.shape[-1]))
    
    for iv in [0,1,3,4,5]:
        if iv == 3:
            sdphi=np.zeros((*RegPolarX.shape,5))
            phi = PARROT[:,:,iv]
            phi[phi<0] = np.nan
            phi[phi>360] = np.nan
            sdphi[4:-2,:,:]=sliding_window_view(phi[2:,:], 5, axis=0)
            interpolator.values = np.nanstd(sdphi,axis = 2, ddof=1).reshape(-1,1)
        else:
            interpolator.values = PARROT[:,:,iv].reshape(-1,1)
        PARITP[:,:,iv] = interpolator(Cx, Cy)
    # scipy.io.savemat('../mat/pyPARROT.mat', {"PARITP": PARITP})

    V_window = sliding_window_view(PARITP[:,:,1], (3, 3))
    V_window = V_window.reshape((*V_window.shape[:2],-1))
    cbr = np.sum(~np.isnan(V_window),axis = 2)/9
    SD_buf = np.zeros(V_window.shape[:2])
    SD_buf[cbr>=0.3] = np.nanstd(V_window[cbr>=0.3].reshape(-1,9),axis = 1, ddof=1)
    stda = np.zeros(Cx.shape)
    stda[1:-1,1:-1] = SD_buf
    
    ########## FTC beta #############
    ############# Beta Cell ############
    a2 = PARITP[:,:,0]
    center_indices = np.argwhere(a2>cellthresh)
    c_indices = clean_indices(center_indices, a2.shape, cellINT)
    cidx = (c_indices[np.newaxis,:,:] + Celldp).astype(int)
    cbox = a2[cidx[:,:,0],cidx[:,:,1]]
    cbr = np.sum( cbox>cellthresh,0)/Celldp.shape[0]
    cbox = cbox[:,cbr>cbcellthrsh]
    c_indices = c_indices[cbr>cbcellthrsh,:]
    llscore = np.zeros(cbox.shape)
    llscore[cbox<=s2xnum[0]] = s2ynum[0]
    pp = np.logical_and(cbox>=s2xnum[0], cbox<s2xnum[1])
    llscore[pp] = s2g*cbox[pp]+s2gc
    llscore[cbox>=s2xnum[1]] = s2ynum[1]
    clscore = np.nansum(llscore,0)/Celldp.shape[0]
    totscore = np.zeros(Cx.shape)
    totscore[c_indices[:,0],c_indices[:,1]] = clscore
    CELLline = medfilt2d(totscore, kernel_size=11)
    
    a2 = CELLline
    center_indices = np.argwhere(a2>cellcsrthresh)
    c_indices = clean_indices(center_indices, a2.shape, widecellINT)
    cidx = (c_indices[np.newaxis,:,:] + Celldp).astype(int)
    cbox = a2[cidx[:,:,0],cidx[:,:,1]]>cellcsrthresh
    cbr = np.sum( cbox>cellcsrthresh,0)/Celldp.shape[0]
    center_indices = c_indices[cbr<1,:]
    cidx = (center_indices[np.newaxis,:,:] + Celldpw).astype(int)
    # a2[cidx[:,:,0],cidx[:,:,1]] = 1
    widecellz = a2>0.5
    ############# Beta Cell ############
    ############# Beta Z, dZ ############
    rotgz = PARITP[:,:,0]
    interpolator.values = diffz.reshape(-1,1)
    rotitp = interpolator(Cx, Cy)
    
    # # cnum1, cnum2, csig1, cfactor1, cintersec1, csig2, cfactor2, cintersec2, cyfill = c_para
    # z_cfun = make_ftc_cscore([15, 20, 3, 3, -1, 12, 4, -2, 3])
    # z_sfun = make_ftc_sscore([0, 5, 5, 2,-1, 5, 3,-3, 1])
    z_cfun = make_ftc_cscore([3, 8, 3, 3, -1, 12, 4, -2, 3])
    z_sfun = make_ftc_sscore([-8, -3, 5, 2,-1, 5, 3,-3, 1])
    zftcs = [FTC_PLAN(ftcc,z_cfun,3), FTC_PLAN(ftcs,z_sfun,1)]
    zbeta = gen_beta(rotgz,thrREF,zftcs)
    # dz_cfun = make_ftc_cscore([5,15,4,3,-2,9,4,-3,2])
    # dz_sfun = make_ftc_sscore([-10,5,5,2,-1,8,2,-3,1])
    dz_cfun = make_ftc_cscore([0,10,4,3,-2,9,4,-3,2])
    dz_sfun = make_ftc_sscore([-10,-5,5,2,-1,8,2,-3,1])
    dzftcs = [FTC_PLAN(ftcc,dz_cfun,2), FTC_PLAN(ftcs,dz_sfun,1)]
    dzbeta = gen_beta(rotitp,thrdREF,dzftcs)
    ############# Beta Z, dZ ############
    zbeta[zbeta<0]=0
    dzbeta[dzbeta<0]=0
    pbeta = (zbeta+dzbeta)/2
    pbeta[np.isnan(PARITP[:,:,0])] = np.nan
    beta = pbeta-widecellz
    beta[beta<0] = 0
    ########## FTC beta #############

    ########## NFGDA eval ###########
    inputNF = np.zeros((*Cx.shape,6))
    inputNF[:,:,0] = beta
    inputNF[:,:,1] = PARITP[:,:,0] # reflectivity
    inputNF[:,:,2] = PARITP[:,:,4] # cross_correlation_ratio
    inputNF[:,:,3] = PARITP[:,:,5] # differential_reflectivity
    inputNF[:,:,4] = stda
    inputNF[:,:,5] = PARITP[:,:,3]

    pnan = np.isnan(inputNF)
    pnansum = np.max(pnan,2)
    inputNF[pnansum,:] = np.nan

    outputGST = fuzzGST.eval_fis(inputNF)
    ########## NFGDA raw output ###########
    ########## post-processing  ###########
    # hh = outputGST>=0.24
    hh = outputGST>=0.6
    hGST = medfilt2d(hh.astype(float), kernel_size=3)
    # smoothedhGST = gaussian_filter(hGST, sigma=1, mode='nearest')
    # skel_nfout = skeletonize(smoothedhGST > 0.3)

    binary_mask = post_moving_avg(hGST) >= 0.6  # Thresholding
    pskel_nfout = binary_dilation(binary_mask, disk(5))
    skel_nfout = skeletonize(pskel_nfout*hh)
    skel_nfout2 = remove_small_objects(skel_nfout, min_size=10, connectivity=2)

    # matout = os.path.join(exp_preds_event,'nf_pred'+os.path.basename(ifn)[5:-3]+'mat')
    data_dict = {"nfout": skel_nfout2,"inputNF":inputNF,
                "timestamp":PARROT_buf["timestamp"]}
    if evalbox_on:
        label_path = os.path.join('../V06/',process_tag,process_tag+'_labels')
        mhandpick = os.path.join(label_path,ifn.split('/')[-1][9:-4]+'.mat')
        try:
            handpick = scipy.io.loadmat(mhandpick)
            evalbox = handpick['evalbox']
        except:
            tprint(ng_tag+f'Warning: No {mhandpick} filling zeros.')
            evalbox = np.zeros(Cx.shape)
        interpolator.values = diffz.reshape(-1,1)
        diffz = interpolator(Cx, Cy)
        data_dict.update({"evalbox":evalbox, \
            "diffz": diffz, \
            'outputGST':outputGST})

    # scipy.io.savemat(matout, data_dict)
    np.savez(nf_path.get_nf_detection_name(ifn,path_config), **data_dict)
    if use_gdal:
        path_config.gdal_writer.log_geo_tif(nf_path.get_nf_detection_name(ifn,path_config,ext='tif'),inputNF)
    nfgda_fig(ifn)

END_GATE = 400
NUM_AZ = 720
var_2_parrot_idx = {'reflectivity': 0, 'velocity': 1, 'spectrum_width': 2, 'differential_phase': 3,
                'cross_correlation_ratio': 4, 'differential_reflectivity': 5}

def get_nexrad(path_config,buf):
    l2_file = buf.filename
    if not os.path.exists(os.path.join(path_config.V06_dir,l2_file)):
        aws_int.download(buf, path_config.V06_dir)
        tprint(dl_tag+f"Got Volume: {l2_file}")
        convert_v06_to_nf_input(l2_file,path_config)
    else:
        tprint(dl_tag+f"Already downloaded. Skip: {l2_file}")
    return

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

def convert_v06_to_nf_input(l2_file, path_config,debug=False):
    v06_file = os.path.join(path_config.V06_dir,l2_file)
    if not (l2_file.endswith('_V06') or l2_file.startswith('._')):
        tprint(cv_tag+"Skip: ", l2_file)
        return

    tprint(cv_tag+"Processing ", l2_file)

    py_path = nf_path.get_nf_input_name(l2_file, path_config)
    # read l2 data
    radar_obj = pyart.io.read_nexrad_archive(v06_file)

    # TODO(pjatau) erase below.
    nsweeps = radar_obj.nsweeps
    vcp = radar_obj.metadata['vcp_pattern']

    debug and tprint(cv_tag+"VCP: ", vcp)

    # VCP 212.
    # slices 0-2-4 contain only dual-pol. super res.
    # slices 1-3-5 contain vel products. super res.
    # slices >= 6 contain all products. normal res.

    # Initialize data cube
    PARROT = np.ma.array(np.ma.array(np.full((END_GATE, NUM_AZ, 6), np.nan, dtype=np.float64)), mask=np.full((END_GATE, NUM_AZ, 6), True))
    in_parrot = np.full(6,False)

    for slice_idx in range(nsweeps):
        radar_range, az_sweep_deg, radar_el, data_slice, mask_slice, labels_slice, data_mask_slice = ReadRadarSliceUpdate(
            radar_obj, slice_idx)
        debug and tprint(cv_tag + "Processing elevation {} degrees".format(np.nanmedian(radar_el)))
        scan_el = np.nanmedian(radar_el)
        # if abs(scan_el-target_el)>0.3:
        #     continue

        i_zero_az = np.argmin(np.abs(az_sweep_deg))
        az_shift = -i_zero_az

        var_idx_slice = {labels_slice[i]: i for i in range(len(labels_slice))}

        for var in var_2_parrot_idx.keys():
            i_var = var_idx_slice[var]
            i_parrot = var_2_parrot_idx[var]

            if not mask_slice[i_var] or in_parrot[i_parrot]:
                continue
            in_parrot[i_parrot] = True

            debug and tprint(cv_tag + "Processing {}. parrot idx {}".format(var, i_parrot))

            curr_data = data_slice[i_var][:, :END_GATE]
            curr_mask = data_mask_slice[i_var][:, :END_GATE]
            # curr_data[curr_mask] = np.nan  # (720, 400)
            curr_data = np.roll(a=curr_data, shift=az_shift, axis=0)
            PARROT[:, :, i_parrot] = curr_data.T
            PARROT[:, :, i_parrot].mask = curr_mask.T
        if np.min(in_parrot):
            timestamp=np.datetime64(pyart.graph.common.generate_radar_time_sweep(radar_obj,slice_idx))
            debug and tprint(cv_tag+"slice idx {} timestamp".format(slice_idx),timestamp)
            break
    # scipy.io.savemat(output_path, {"PARROT": PARROT})
    np.savez(py_path, PARROT=PARROT.data,mask=PARROT.mask,timestamp=timestamp)

GFG_NONE    = 0b00000
GFG_LABEL   = 0b00001
GFG_NFGDA   = 0b00010
GFG_PREDICT = 0b00100
class GFGroups:
    ogn = np.array([Cx[0,0],Cy[0,0]])[:,np.newaxis]
    shape = Cx.shape
    n_anchors = 10
    dummy_anchors = np.array([np.linspace(Cx[0,0],Cx[0,1],10), np.linspace(Cy[0,0],Cy[0,1],10)])[np.newaxis,:]
    def __init__(self, arc_anchors = None, timestamp=None, datakind = GFG_NONE):
        if arc_anchors is None:
            self.arc_anchors = self.dummy_anchors.copy()
        else:
            self.arc_anchors = np.asarray(arc_anchors).copy()
        self.timestamp = timestamp
        self.cur_motions = np.full(self.arc_anchors.shape, np.nan)
        self.pre_motions = np.full(self.arc_anchors.shape, np.nan)
        self.next_gp = np.full(self.arc_anchors.shape[0], np.nan)
        self.datakind = datakind
    def anchors_to_arcs(self):
        self.arc_points = []
        if len(self.arc_anchors)!=0:
            for nf_anchor in self.arc_anchors:

                self.arc_points.append(rotation_polyfit(nf_anchor,2))
        return self.arc_points
    def anchors_to_arcs_map(self):
        if len(self.arc_anchors)!=0:
            self.anchors_to_arcs()
            fit_points = np.concatenate(self.arc_points, axis=1)
            self.arcs_map = points_to_binary_grid(fit_points,self.shape,self.ogn,0.5)
            return self.arcs_map
        else:
            return np.zeros(self.shape,dtype=bool)

    def save(self, file_path: str | Path):
        """Save the GFGroups data to an .npz file"""
        np.savez(
            file_path,
            arc_anchors=self.arc_anchors,
            timestamp=self.timestamp,
            cur_motions=self.cur_motions,
            pre_motions=self.pre_motions,
            next_gp=self.next_gp,
            datakind = self.datakind,
        )

    @classmethod
    def load(cls, file_path: str | Path) -> "GFGroups":
        """Load a GFGroups object from an .npz file"""
        file_path = Path(file_path)
        data = np.load(file_path, allow_pickle=True)
        obj = cls.__new__(cls)  # create instance without calling __init__
        # assign directly
        obj.arc_anchors = data["arc_anchors"]
        obj.timestamp = data["timestamp"]
        obj.cur_motions = data["cur_motions"]
        obj.pre_motions = data["pre_motions"]
        obj.next_gp = data["next_gp"]
        obj.datakind = data["datakind"]
        return obj

class DataGFG(GFGroups):
    def __init__(self, data, binary_mask, kind = GFG_NONE):
        arc_anchors = []
        neighbors = convolve(binary_mask.astype(int), kernel, mode='constant') - binary_mask
        branch = (binary_mask & (neighbors >= 3))
        # Temporarily remove branch points
        skel_wo_branches = binary_mask.copy()
        skel_wo_branches[branch] = 0
        self.groups = label(skel_wo_branches, connectivity=2)
        # self.groups = label(binary_mask, connectivity=2)
        reduce = []
        for im in range(1,np.max(self.groups)+1):
            mask = self.groups == im
            gf_points = np.array([Cx[mask],Cy[mask]])
            if gf_points.shape[1]>2:
                rot_points = rotation_polyfit(gf_points,2,self.n_anchors)
                arc_anchors.append(rot_points)
            else:
                reduce.append(im)
        for im in reduce[::-1]:
            self.groups[self.groups==im] = 0
            self.groups[self.groups>=im] -= 1
        super().__init__(np.array(arc_anchors),data['timestamp'],kind)

class Prediction_Connection:
    def __init__(self, endidx, flip, igps, egps):
        self.endidx = endidx
        self.flip = flip
        self.igp_anchor = np.array(igps.arc_anchors)
        self.egp_anchor = np.array(egps.arc_anchors)
        self.motions = np.full(self.igp_anchor.shape,np.nan)
        self.speeds = np.full(self.igp_anchor.shape[0],np.nan)
        self.dt = (egps.timestamp - igps.timestamp)/np.timedelta64(60, 's')
        self.make_motion()
        self.copy_motion(igps, egps)
    def make_motion(self):
        for ii, ie in enumerate(self.endidx):
            ep = self.egp_anchor[ie]
            if self.flip[ii]:
                ep = np.fliplr(ep)
            self.motions[ii] = ep-self.igp_anchor[ii]
        self.motions = self.motions/self.dt 
        if self.motions.ndim>1:
            self.speeds = np.sqrt(np.sum(self.motions**2,axis=1))
    def copy_motion(self, igps, egps):
        igps.next_gp = self.endidx
        igps.cur_motions = self.motions
        for ii, ie in enumerate(set(self.endidx)):
            egps.pre_motions[ie] = np.mean(self.motions[self.endidx==ie].reshape(-1,2,egps.n_anchors),axis=0)
            # if np.sum(self.endidx==ie)>1:
            #     print(egps.pre_motions[ie],self.motions[self.endidx==ie].shape,self.flip[self.endidx==ie])

class Prediction_Worker:
    def __init__(self, gps):
        self.gps = gps
        self.connects = {}
    def update_velocitys(self, k=None):
        if k is None:
            for ig in range(len(self.gps)-1):
                A = self.gps[ig]
                B = self.gps[ig+1]
                self.GFG_motion(A,B,ig)
        else:
            A = self.gps[k]
            B = self.gps[k+1]
            self.GFG_motion(A,B,k)
    def GFG_motion(self,sgfg,egfg,conn):
        A = sgfg.arc_anchors
        B = egfg.arc_anchors
        if (A.ndim != 3) or (B.ndim != 3):
            # print(A,B)
            self.connects[conn] = Prediction_Connection([], [], sgfg, egfg)
        else:
            A_exp = A[:, None, :, :]
            B_exp = B[None, :, :, :]

            # Compute normal distances (mean over last axis=2, then mean over coordinates)
            dist_norm = np.linalg.norm(A_exp - B_exp, axis=2)
            
            # Compute flipped distances
            B_flip = B[:, :, ::-1]  # flip along the sample axis (10)
            B_flip_exp = B_flip[None, :, :, :]
            dist_flip = np.linalg.norm(A_exp - B_flip_exp, axis=2)
            
            # Choose the smaller distance
            dist_final = np.max(dist_norm, axis=-1).copy()
            flip_arc = np.sum(dist_norm, axis=-1) > np.sum(dist_flip, axis=-1)
            dist_final[flip_arc]=np.max(dist_flip, axis=-1)[flip_arc]
            
            endpoint = np.argmin(dist_final,1)
            self.connects[conn] = Prediction_Connection(endpoint, flip_arc[np.arange(endpoint.size),endpoint], sgfg, egfg)
            # return endpoint, flip_arc[np.arange(endpoint.size),endpoint], dist_final[np.arange(endpoint.size),endpoint]
    def prediction(self,startf,dt,mode='element'):
        if mode == 'mean':
            anchors = self.connects[startf].igp_anchor+np.mean(self.connects[startf].motions,axis=2)[:,:,np.newaxis]*dt
            tstp = self.gps[startf].timestamp + dt*np.timedelta64(1, 'm')
        elif mode == 'element':
            anchors = self.connects[startf].igp_anchor+self.connects[startf].motions*dt
            tstp = self.gps[startf].timestamp + dt*np.timedelta64(1, 'm')
        return GFGroups(anchors,tstp)

def eval_nf(nfloc,evalline,evalbox,gd):
    nfpredict = dilation(nfloc, disk(5))
    Mhits = np.logical_and(evalline,nfpredict)
    Mmiss = np.logical_and(evalline,~nfpredict)
    q_arc = nf_arc(gd['Cx'],gd['Cy'],nfloc)
    Mfa = np.logical_and(np.logical_not(evalbox),q_arc)
    sch,scm,scf,scp = np.sum(Mhits), np.sum(Mmiss), np.sum(Mfa), np.sum(q_arc)
    HR = 1e2*sch/(scm+sch)
    FAR = 1e2*scf/(scp)
    return (HR, FAR), (Mhits, Mmiss, Mfa, q_arc)

def log_stat(fn, stats_list):
    grouped = list(zip(*stats_list))
    # Stack each group
    stacked = [np.stack(arr_list) for arr_list in grouped]
    # Create a dictionary with key names
    save_dict = {
    'hits':stacked[0], 
    'miss':stacked[1], 
    'fa':stacked[2], 
    'q_arc':stacked[3]
    }
    # Save to .npz
    np.savez(fn, **save_dict)
    return save_dict

###### Polyfit Arc Functions ##########
def rotation_matrix_2d(theta):
    """
    Returns a 2×2 rotation matrix for rotating points counterclockwise by angle theta (in radians).
    """
    return np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])

def find_roation_coord(points):
    points_T = points.T  # shape = (N, 2)
    dist_matrix = squareform(pdist(points_T))
    i, j = np.unravel_index(np.argmax(dist_matrix), dist_matrix.shape)
    max_dist = dist_matrix[i, j]
    delta = points[:,i] - points[:,j]
    angle = np.arctan2(delta[1], delta[0])
    return angle, points[:,i]

def rotation_polyfit(points,n,fitn=None):
    if np.isnan(points).any():
        return points
    if points.shape[1]>2:
        try:
            ang, origin = find_roation_coord(points)
            rot_points = np.matmul(rotation_matrix_2d(-ang),points-origin[:,np.newaxis])
            coeffs = np.polyfit(rot_points[0,:], rot_points[1,:], n)
        except:
            return points
        if fitn is None:
            fx = np.arange(np.min(rot_points[0,:]),np.max(rot_points[0,:])+0.25,0.25)
        else:
            fx = np.linspace(np.min(rot_points[0,:]),np.max(rot_points[0,:]),fitn)
        fy = np.polyval(coeffs, fx)
        return np.matmul(rotation_matrix_2d(ang),np.array([fx,fy]))+origin[:,np.newaxis]
    else:
        return points

def points_to_binary_grid(xy_points, gshape, origin, dg):
    grid = np.zeros(gshape, dtype=bool)

    # Round or floor the points to nearest integers
    xy_int = ((xy_points-origin)//dg).astype(int)

    # Clip to ensure they fall within bounds
    xy_int[0, :] = np.clip(xy_int[0, :], 0, gshape[1] - 1)
    xy_int[1, :] = np.clip(xy_int[1, :], 0, gshape[0] - 1)

    grid[xy_int[1, :], xy_int[0, :]] = True # grid [y, x]
    return grid

def nf_arc(xx,yy,bw):
    groups = label(bw, connectivity=2)
    fit_chunks = []
    for im in range(1,np.max(groups)+1):
        mask = groups == im
        gf_points = np.array([xx[mask],yy[mask]])
        rot_points = rotation_polyfit(gf_points,2)
        fit_chunks.append(rot_points)
    ogn = np.array([xx[0,0],yy[0,0]])[:,np.newaxis]
    if len(fit_chunks)!=0:
        fit_points = np.concatenate(fit_chunks, axis=1)
        return points_to_binary_grid(fit_points,xx.shape,ogn,0.5)
    else:
        return np.zeros(xx.shape,dtype=bool)
###### Polyfit Arc Functions ##########

# def post_proc(inGST):
#     hGST = medfilt2d(inGST.astype(float), kernel_size=3)
#     binary_mask = post_moving_avg(hGST) >= 0.6  # Thresholding
#     pskel_nfout = binary_dilation(binary_mask, disk(5))
#     skel_nfout = skeletonize(pskel_nfout*inGST)
#     skel_nfout2 = remove_small_objects(skel_nfout, min_size=10, connectivity=2)
#     return skel_nfout2

# gfv = [4, 32]
# class GFSpace:
#     def __init__(self, gfv = [18, 72]):
#         self.gfv = gfv
#         self.data =[]
#         self.nf_pair = []
#         self.nf_loc = []
#         self.tstamp = []
    
#     def load_nf(self,fn):
#         self.data.append(np.load(fn))
#         self.tstamp.append(get_tstamp(fn))
#         self.nf_loc.append( {'x':Cx[self.data[-1]['nfout']],'y':Cy[self.data[-1]['nfout']],
#                              'idx':np.arange(Cx.size).reshape(Cx.shape)[self.data[-1]['nfout']]})
#         if len(self.nf_loc)>1:
#             self.nf_pair.append(self.connect_nf(self.nf_loc[-2],self.nf_loc[-1],(self.tstamp[-1]-self.tstamp[-2]).total_seconds()))
#         else:
#             self.shp = Cx.shape
            
#     def connect_nf(self,t1,t2,dt):
#         A = np.column_stack((t1['x'].reshape(-1), t1['y'].reshape(-1)))
#         B = np.column_stack((t2['x'].reshape(-1), t2['y'].reshape(-1)))
#         # coord [dim_sample, dim_xy]
#         # coord [dim_sample, 2]
#         tree = cKDTree(A)
#         dists_for, indices_for = tree.query(B)
#         tree = cKDTree(B)
#         dists_bac, indices_bac = tree.query(A)
#         pair_pool = np.vstack((np.column_stack((indices_for,np.arange(B.shape[0]))),np.column_stack((np.arange(A.shape[0]),indices_bac)))).astype(int)
#         dists_pool = np.concatenate((dists_for,dists_bac),axis=0)
#         mask = np.logical_and(dists_pool > self.gfv[0]*dt/3600, dists_pool < self.gfv[1]*dt/3600)
#         return np.concatenate((pair_pool[mask,:],np.zeros(np.sum(mask),dtype=bool).reshape(-1,1)),axis=1)
        
#         # return pair_pool[mask,:].astype(int)
#     def clean_short_track(self):
#         for ic in range(len(self.nf_pair)-1):
#             keep_mask = np.logical_or(np.isin(self.nf_pair[ic][:,1], np.unique(self.nf_pair[ic+1][:,0])),self.nf_pair[ic][:,2]>0)
#             self.nf_pair[ic] = self.nf_pair[ic][keep_mask,:]
#             self.nf_pair[ic+1][:,2] = np.logical_or(self.nf_pair[ic+1][:,2], np.isin(self.nf_pair[ic+1][:, 0], self.nf_pair[ic][:, 1]))
    
#     def clean_random_track_motion(self):
#         self.cal_motion()
#         for ic in range(len(self.nf_pair)-1):
#             curdir = self.nf_pair[ic][:,-2]+1j*self.nf_pair[ic][:,-1]
#             nextdir = np.zeros(curdir.size)
#             for ip in range(curdir.size):
#                 mask = self.nf_pair[ic+1][:,0]==self.nf_pair[ic][ip,1]
#                 nextdir[ip] = np.mean(self.nf_pair[ic+1][mask,-2]+1j*self.nf_pair[ic+1][mask,-1],axis=0)
#             dirdiff = get_dirdiff(curdir,nextdir)
#             self.nf_pair[ic]=np.concatenate((self.nf_pair[ic],dirdiff.reshape(-1,1)),axis=1)
        
#     def cal_motion(self):
#         for ic in range(len(self.nf_pair)):
#             t1 = self.nf_loc[ic]
#             t2 = self.nf_loc[ic+1]
#             A = np.column_stack((t1['x'].reshape(-1), t1['y'].reshape(-1)))
#             B = np.column_stack((t2['x'].reshape(-1), t2['y'].reshape(-1)))
#             start = A[self.nf_pair[ic][:,0]]
#             end = B[self.nf_pair[ic][:,1]]
#             motion = end-start
#             self.nf_pair[ic]=np.concatenate((self.nf_pair[ic],motion),axis=1)
    
#     def get_cln_nf(self,ic):
#         buf = np.zeros(self.shp,dtype=bool).reshape(-1)
#         buf[self.nf_loc[ic]['idx'][self.nf_pair[ic][:,0].astype(int)]]=True
#         return buf.reshape(self.shp)
#         # return remove_small_objects(buf.reshape(self.shp), min_size=5, connectivity=2)

def get_tstamp(ppi_file):
    ppi_id = os.path.basename(ppi_file)
    ppi_name = ppi_id[11:]
    date_part = ppi_name[4:12]   # 5:12 in MATLAB → 4:12 in Python
    time_part = ppi_name[13:19]  # 14:19 in MATLAB → 13:19 in Python
    tstamp_date = datetime.strptime(date_part, "%Y%m%d")
    tstamp_time = datetime.strptime(time_part, "%H%M%S").time()
    tstamp = datetime.combine(tstamp_date.date(), tstamp_time)
    return tstamp

def get_dirdiff(dir1,dir2):
    return np.rad2deg(np.angle(dir1*np.conj(dir2)))

def nfgda_fig(l2_file):
    py_path = nf_path.get_nf_detection_name(l2_file, path_config)
    
    data = np.load(py_path)
    gps = DataGFG(data,data['nfout'])
    fig, axs = plt.subplots(1, 1, figsize=(3.3/0.7, 3/0.7),dpi=250)
    pdata = np.ma.masked_where(rmask,data['inputNF'][:,:,1])
    pcz=axs.pcolormesh(Cx,Cy,pdata,cmap=cl.zmap,norm=cl.znorm)

    # if np.sum(data['nfout'])>0:
    if gps.arc_anchors.ndim==3:
        axs.plot(gps.arc_anchors[:,0,:].T,gps.arc_anchors[:,1,:].T,alpha=0.7,color='k')

    valid_time = gps.timestamp
    fig.suptitle(valid_time.astype(datetime.datetime).item().strftime('%Y/%m/%d %H:%M:%S'),y=0.95)
    axs.set_xlim(-100,100)
    axs.set_ylim(-100,100)
    axs.set_xlabel('x(km)')
    axs.set_ylabel('y(km)',labelpad=-10)
    axs.set_aspect('equal')
    fig.savefig(py_path[:-3]+'png')
    plt.close(fig)

def nfgda_forecast(l2_file_0,l2_file_1,debug=False,suppress_fig=False):
    if suppress_fig:
        tprint(df_tag+
            f'{l2_file_0} figures are suppressed.')
    py_path = nf_path.get_nf_detection_name(l2_file_0, path_config)
    data = np.load(py_path)
    py_path = nf_path.get_nf_detection_name(l2_file_1, path_config)
    data1 = np.load(py_path)
    gps=[DataGFG(data,data['nfout']),DataGFG(data1,data1['nfout'])]
    worker = Prediction_Worker(gps)
    worker.update_velocitys(0)
    tvec = worker.gps[0].timestamp.astype('datetime64[m]')\
        +np.arange(60,7201,60)*np.timedelta64(1, 's')
    fig, axs = plt.subplots(1, 1, figsize=(3.3/0.7, 3/0.7),dpi=150)
    pdata = np.ma.masked_where(rmask,data['inputNF'][:,:,1])
    pcz=axs.pcolormesh(Cx,Cy,pdata,cmap=cl.zmap,norm=cl.znorm)
    axs.set_xlim(-100,100)
    axs.set_ylim(-100,100)
    axs.set_xlabel('x(km)')
    axs.set_ylabel('y(km)',labelpad=-10)
    axs.set_aspect('equal')
    forecast_anchors = []
    for t in tvec:
        dt = (t-worker.gps[0].timestamp)/np.timedelta64(60, 's')
        if worker.connects[0].motions.ndim==3:
            ende = worker.prediction(0, dt)
            endm = worker.prediction(0, dt,mode='mean')
            forecast_anchors.append((ende,endm))
            if suppress_fig: continue
            axs.plot(ende.arc_anchors[:,0,:].T,ende.arc_anchors[:,1,:].T,alpha=0.7,color='k')
            axs.plot(endm.arc_anchors[:,0,:].T,endm.arc_anchors[:,1,:].T,alpha=0.7,color='r')
        else:
            forecast_anchors.append((GFGroups(timestamp=t),GFGroups(timestamp=t)))
            debug and tprint(df_tag+f'{C.RED_B} Prediction dimension != 3 {C.RESET}',worker.connects[0].motions)
        if suppress_fig: continue
        fig.suptitle(worker.gps[0].timestamp.astype(datetime.datetime).item().strftime('%Y/%m/%d %H:%M:%S')+'\n'+
                t.astype(datetime.datetime).strftime('%Y/%m/%d %H:%M:%S')+f' (+{int(dt)} mins)',y=0.97)
        fig.savefig(nf_path.get_nf_forecast_name(l2_file_0, path_config,t))
        for ln in axs.lines[:]:
            ln.remove()
    plt.close(fig)
    with open(nf_path.get_nf_forecast_pkl_name(l2_file_0, path_config), "wb") as f:
        pickle.dump((tvec,forecast_anchors), f)
    return tvec,forecast_anchors

    # ele_map=[]
    # mean_map=[]
    # for t in tvec:
    #     dt = (t-worker.gps[0].timestamp)/np.timedelta64(60, 's')
    #     end = worker.prediction(0, dt)
    #     ele_map.append(end.anchors_to_arcs_map())
    #     end = worker.prediction(0, dt,mode='mean')
    #     axs.plot(gps.arc_anchors[:,0,:].T,gps.arc_anchors[:,1,:].T,alpha=0.7,color='k')
    #     mean_map.append(end.anchors_to_arcs_map())
    # data_dict = {"ele_map": ele_map, "mean_map": mean_map, "start_timestamps":tvec}
    # np.savez(nf_path.get_nf_forecast_name(l2_file_0, path_config), **data_dict)

def exp_weight(dt,t_const):
    return np.exp(-dt/t_const)

def nfgda_stochastic_summary(forecasts,l2_file_0,force=False):
    py_path = nf_path.get_nf_detection_name(l2_file_0, path_config)
    data = np.load(py_path)

    live_tdx = np.full((len(forecasts),),False)
    tstart = data['timestamp']
    tnow = tstart.astype('datetime64[m]')+60*np.timedelta64(1, 's')
    forecast_end = tstart
    for tdx,buf in enumerate(forecasts):
        if len(buf)==2:
            if buf[0][0]>tnow and not(force):
                tprint(sf_tag+
                    f'{tnow} is out of date. New data {buf[0][0]} available.')
                return
            if buf[0][-1]>tnow and buf[0][0]<tnow:
                live_tdx[tdx]=True
            forecast_size = buf[0].size
            # tprint(sf_tag+f'forecast_size = {forecast_size}')
    if np.sum(live_tdx)==0:
        tprint(sf_tag+
            'No forecast for summary.')
        return
    else:
        summary_tdx = np.where(live_tdx)[0]
        tprint(sf_tag+
            f'Summary forecasts[{len(summary_tdx)}]:',
            *[forecasts[tdx][0][0] for tdx in summary_tdx])
    ips=[]
    for tdx in summary_tdx:
        ips.append(np.where(forecasts[tdx][0]==tnow)[0][0])
    ips = np.array(ips,dtype=int)
    tline = forecasts[summary_tdx[np.argmin(ips)]][0][np.min(ips):]
    fig, axs = plt.subplots(1, 1, figsize=(3.3/0.7, 3/0.7),dpi=150)
    pdata = np.ma.masked_where(rmask,data['inputNF'][:,:,1])
    pcz=axs.pcolormesh(Cx,Cy,pdata,cmap=cl.zmap,norm=cl.znorm)
    axs.set_xlim(-100,100)
    axs.set_ylim(-100,100)
    axs.set_xlabel('x(km)')
    axs.set_ylabel('y(km)',labelpad=-10)
    axs.set_aspect('equal')
    # for ipdx,tdx in zip(ips,summary_tdx):
    #     print(forecasts[tdx][0][ipdx])

    for ip in range(forecast_size-np.min(ips)):
        ps = 0
        pgf = np.zeros(Cx.shape)
        valid_time = tline[ip]
        for ipdx,tdx in zip(ips,summary_tdx):
            if (ipdx+ip) >= forecast_size:continue
            arcs = forecasts[tdx][1][ipdx+ip]
            if valid_time!= forecasts[tdx][0][ipdx+ip]:
                raise ValueError(sf_tag+f'sf time mismatch {valid_time} {forecasts[tdx][0][ipdx+ip]}')
            dt = (valid_time-tnow)/np.timedelta64(60, 's')
            ele_w = exp_weight(dt,ele_t_const)
            mean_w = exp_weight(dt,mean_t_const)
            ps += ele_w + mean_w
            pgf += ele_w*binary_dilation(arcs[0].anchors_to_arcs_map(), footprint=disk(3)).astype(float)
            pgf += mean_w*binary_dilation(arcs[1].anchors_to_arcs_map(), footprint=disk(3)).astype(float)
        pgf=pgf/ps*1e2
        pgf[rmask] = 0
        
        fig.suptitle(tnow.astype(datetime.datetime).strftime('%Y/%m/%d %H:%M:%S')+'\n'
            +valid_time.astype(datetime.datetime).strftime('%Y/%m/%d %H:%M:%S')+f' (+{int(dt)} mins)',y=0.97)
        cs = axs.contour(Cx, Cy, pgf, levels=[30],colors='red')
        fig.savefig(nf_path.get_nf_s_forecast_name(path_config,valid_time))
        cs.remove()
        data_dict = {"nfproxy": pgf, "timestamp":valid_time}
        np.savez(nf_path.get_nf_s_forecast_name(path_config,valid_time)[:-3]+'npz', **data_dict)
