import numpy as np
from matplotlib import pyplot as plt
'''import yt
import pyxsim
import soxs
import moviepy
import moviepy.video.io.ImageSequenceClip
from moviepy.editor import *'''

import sys
sys.path.append('/home/nn0933/athenaresearch/python/') # path to athena_data and athena_kit module
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["figure.facecolor"] = "white"
from athena_data_2023_01 import AthenaBinary
from athena_data_2023_01 import AthenaBinaries
from athena_data_2023_01 import load_dict_from_hdf5, save_dict_to_hdf5
import athena_kit as ak

mu=0.618
muH=1.4
unit=ak.Units(lunit=ak.pc_cgs,munit=mu*ak.atomic_mass_unit_cgs*ak.pc_cgs**3,mu=mu)
path='/home/mg9443/scratch/share/data/snr/snr_net/snrnet_h05_t16_128_k0/data/bin/'
atbs=AthenaBinaries(binarypath=path+f"/Blast.hydro_w.",path=path,version='0.2')
files= sorted(glob.glob(path + "*Blast.hydro_w.*"))
fstart= int(files[0][-6:-4])
fend= int(files[-1][-6:-4])
ilist = list(range(fstart, fend + 1))
atbs.read(ilist) # read binaries from these specific files, check number of files 
atbs.config()

x1min=float(abin.header('mesh','x1min'))
x2min=float(abin.header('mesh','x2min'))
x3min=float(abin.header('mesh','x3min'))
x1max=float(abin.header('mesh','x1max'))
x2max=float(abin.header('mesh','x2max'))
x3max=float(abin.header('mesh','x3max'))
xyz=np.array([x1min,x1max,x2min,x2max,x3min,x3max])
dens_arr=abin.get_data('dens',xyz=list(xyz*1),level=0).T #to account for Athena's zyx
temp_arr=abin.get_data('temp',xyz=list(xyz*1),level=0).T
velx_arr=abin.get_data('velx',xyz=list(xyz*1),level=0).T
vely_arr=abin.get_data('vely',xyz=list(xyz*1),level=0).T
velz_arr=abin.get_data('velz',xyz=list(xyz*1),level=0).T
pres_arr=abin.get_data('pres',xyz=list(xyz*1),level=0).T
data = dict(density = (dens_arr*unit.density_cgs, "g*cm**-3"), 
            number_density = (dens_arr, "cm**-3"), 
            H_nuclei_density = (mu/muH*dens_arr, "cm**-3"), ##check with Chang-Goo
            temperature= (temp_arr*unit.temperature_cgs, "K"),
            pressure= (pres_arr*unit.pressure_cgs, "erg*cm**-3"),            
            velocity_x = (velx_arr*unit.velocity_cgs, "cm*s**-1"),
            velocity_y = (vely_arr*unit.velocity_cgs, "cm*s**-1"),
            velocity_z = (velz_arr*unit.velocity_cgs, "cm*s**-1"),
           )
bbox = np.array([[x1min,x1max], [x2min,x2max], [x3min,x3max]])
ds = yt.load_uniform_grid(data, dens_arr.shape, length_unit="pc", periodicity=(False,False,False),
                          bbox=bbox, nprocs=256,default_species_fields='ionized')
#Add Xray field
source_model = pyxsim.CIESourceModel("spex", 0.05, 11.0, 1000, 1.0, binscale="log")
xray_fields = source_model.make_source_fields(ds, 0.5, 7.0) 
slc = ds.slice('z',0)
sfrb = slc.to_frb((64,'pc'), (ds.domain_dimensions[0],ds.domain_dimensions[0]))
trial= sfrb[('gas','density')]
slice = dict()
for k in ['density','temperature']:
    slice[k]=sfrb[('gas',k)]
save_dict_to_hdf5(slice,'tmp.hdf5')