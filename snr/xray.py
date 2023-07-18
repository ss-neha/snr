import yt
import numpy as np
import pyxsim
#import soxs
import matplotlib.pyplot as plt
import glob
import sys
import os
import pickle
sys.path.append('/home/nn0933/athenaresearch/python/') # path to athena_data and athena_kit module
sys.path.append('/home/nn0933/athenaresearch/') # path to athena_data and athena_kit module
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["figure.facecolor"] = "white"
from athena_data_2023_01 import AthenaBinary, save_dict_to_hdf5
# from athena_data_2023_01 import AthenaBinaries
# from athena_data_2023_01 import load_dict_from_hdf5, 
import athena_kit as ak

# define heating rate temporarily
hrate = 2e-26

# yt unit class
u = yt.units

# defining extra field functions
def _pok(field, data):
    return (data[("gas", "H_nuclei_density")]
        * data[("gas", "temperature")]
    )

def cooling_rate(field, data):
    nh = ('gas','number_density')
    T = ('gas','temperature')
    nH = data[nh].v
    return nH*nH*ak.CoolFnShure(data[T].v)*u.erg/u.cm**3/u.s   #.unit("erg/cm**(-3)*s**(-1)")

def net_cooling_rate(field, data):
    nh = ('gas','number_density')
    T = ('gas','temperature')
    nH=data[nh].v
    return nH*(nH*ak.CoolFnShure(data[T].v)-hrate)*u.erg/u.cm**3/u.s   #.unt("erg/cm**(-3)*s**(-1)")

def e_density(field, data):
    gamma=5/3.
    return data[("gas", "pressure")]/(gamma-1)


def total_e_density(field, data):
    return (data[("gas", "thermal_energy_density")] +
            data[("gas", "kinetic_energy_density")])

def total_e(field, data):
    return (data[("gas", "total_energy_density")] *
        data[("gas", "cell_volume")])

def energy_loss(field,data):
    return (data[("gas", "net_cooling_rate")] *
            data[("gas", "cell_volume")])

# adding new fields
def add_yt_fields(ds):    
    ds.add_field(
        name=("gas", "pok"),
        function=_pok,
        sampling_type="local",
    )

    ds.add_field(
        name=("gas", "net_cooling_rate"),
        function=net_cooling_rate,
        sampling_type="local",
    )

    ds.add_field(
        name=("gas", "cooling_rate"),
        function=cooling_rate,
        sampling_type="local",
    )

    ds.add_field(
        name=("gas", "thermal_energy_density"),
        function=e_density,
        sampling_type="local",
    )

    ds.add_field(
        name=("gas", "total_energy_density"),
        function=total_e_density,
        sampling_type="local",
    )

    ds.add_field(
        name=("gas", "total_energy"),
        function=total_e,
        sampling_type="local",
    )

    ds.add_field(
        name=("gas", "energy_loss_rate"),
        function=energy_loss,
        sampling_type="local",
    )

def get_edges(bin_info):
    bmin,bmax,N = bin_info['min'],bin_info['max'],bin_info['Nedge']
    if bin_info['log']:
        edges = np.logspace(bmin,bmax,N)
    else:
        edges = np.linspace(bmin,bmax,N)
    return edges

# creating/saving one joint PDF
def make_one_PDF(ds,xf,yf,bin_fields,
                 wflist=['Lx','vol','mass','etot','eloss'],
                 inum=0,
                 save=False, outdir='./final_test/'):
    pdf = dict()
    for wf in wflist:
        xdata = ds.r[bin_fields[xf]['fieldname']].to(bin_fields[xf]['units'])
        ydata = ds.r[bin_fields[yf]['fieldname']].to(bin_fields[yf]['units'])
        wdata = ds.r[bin_fields[wf]['fieldname']].to(bin_fields[wf]['units'])
        xedges = get_edges(bin_fields[xf])
        yedges = get_edges(bin_fields[yf])
        h2d,_,_ = np.histogram2d(xdata, ydata,
                                bins=[xedges,yedges],weights=wdata)
        pdf[wf] = h2d
        pdf[f'{wf}_sum'] = wdata.sum().v
    pdf['x_edges']= xedges
    pdf['y_edges']= yedges
    pdf['xf'] = xf
    pdf['yf'] = yf

    if save:
        fnamebase = f'pdf_{xf}_{yf}_{inum:04d}.p'
        with open(os.path.join(outdir,fnamebase),'wb') as fp:
            pickle.dump(pdf,fp)
    return pdf

def setup_bin_field_info(Nbin=256, xray_field=('gas', 'xray_luminosity_0.5_7.0_keV')):
    # bin field definitions
    bin_fields = dict()
    bin_fields['T'] = dict(fieldname = ('gas','temperature'),
                            units = 'K',
                            Nedge = Nbin+1, min=1, max=9, log=True)
    bin_fields['nH'] = dict(fieldname = ('gas', 'H_nuclei_density'),
                            units = 'cm**-3',
                            Nedge = Nbin+1, min=-3, max=3, log=True)
    bin_fields['pok'] = dict(fieldname = ('gas', 'pok'),
                            units = 'cm**-3*K',
                            Nedge = Nbin+1, min=1, max=8, log=True)
    bin_fields['r'] = dict(fieldname = ('index', 'spherical_radius'),
                    units = 'pc',
                    Nedge = Nbin+1, min=0, max=32, log=False)
    bin_fields['R'] = dict(fieldname = ('index', 'cylindrical_radius'),
                        units = 'pc',
                        Nedge = Nbin+1, min=0, max=32, log=False)
    bin_fields['Lx'] = dict(fieldname = xray_field,
                    units = 'erg/s',
                    Nedge = Nbin+1, min=25, max=40, log=True)
    bin_fields['velz'] = dict(fieldname = ('gas','velocity_z'),
                units = 'km/s',
                Nedge = Nbin+1, min=-1500, max=1500, log=False)
    bin_fields['x'] = dict(fieldname = ('gas','x'),
                units = 'pc',
                Nedge = Nbin+1, min=-32, max=32, log=False)
    bin_fields['y'] = dict(fieldname = ('gas','y'),
            units = 'pc',
            Nedge = Nbin+1, min=-32, max=32, log=False)
            
 
    # weight field only
    bin_fields['vol'] = dict(fieldname = ('gas', 'cell_volume'),
                             units = 'pc**3')
    
    bin_fields['mass'] = dict(fieldname = ('gas', 'cell_mass'),
                             units = 'Msun')    
    bin_fields['etot'] = dict(fieldname = ('gas', 'total_energy'),
                             units = 'erg')    
    bin_fields['eloss'] = dict(fieldname = ('gas', 'energy_loss_rate'),
                             units = 'erg/s')    

    # T=ds.r[('gas','temperature')]
    # vz=ds.r[('gas','velocity_z')]
    # vr=ds.r[('gas','radial_velocity')]
    # vsr=ds.r[('gas', 'velocity_spherical_radius')]
    # vcr=ds.r[('gas', 'velocity_cylindrical_radius')]
    # nH=ds.r[('gas', 'H_nuclei_density')]
    # Lx=ds.r[xray_fields[1]]
    # vol=ds.r[('gas','cell_volume')]
    # mass=ds.r[('gas','cell_mass')]
    # pres=ds.r[('gas','pressure')]
    # pok=ds.r[('gas','pok')]
    # dens=ds.r[('gas','density')]
    # sr=ds.r[('index','spherical_radius')]
    # cr=ds.r[('index','cylindrical_radius')]
    # cool=ds.r[('gas','cooling_rate')]
    # ncool=ds.r[('gas','net_cooling_rate')]
    # etot=ds.r[('gas','total_energy')]
    # eloss=ds.r[('gas','energy_loss_rate')]



    # # bin defnitions
    # T_edges = np.logspace(1,9,Nbin+1)
    # vz_edges = np.linspace(-1500,1500,Nbin+1)
    # Lx=ds.r[xray_fields[1]]
    # mass_edges=np.logspace(28,40,Nbin+1)
    # nH_edges=np.logspace(-3,3,Nbin+1)
    # pres_edges=np.logspace(-17,-8,Nbin+1)
    # vsr_edges=np.linspace(-600,1500,Nbin+1)
    # vcr_edges=np.linspace(-600,1500,Nbin+1)
    # vr_edges=np.linspace(-600,1500,Nbin+1)
    # pok_edges=np.logspace(1,8,Nbin+1)
    # ncool_edges=np.linspace(-5e-24,7e-19,Nbin+1)
    # cr_edges=np.linspace(0,64,Nbin+1)
    # sr_edges=np.linspace(0,64,Nbin+1)
    # Lx_edges=np.logspace(25,35,Nbin+1)
   
    return bin_fields

def ytload(abin):
    # unit initialization
    mu=0.618
    muH=1.4
    unit=ak.Units(lunit=ak.pc_cgs,munit=mu*ak.atomic_mass_unit_cgs*ak.pc_cgs**3,mu=mu)

    # read Athenabinary into yt
    # abin=atbs.abins[inum]
    
    # update global variable hrate
    hrate=float(abin.header("problem","hrate"))

    # local varaibles
    x1min=float(abin.header('mesh','x1min'))
    x2min=float(abin.header('mesh','x2min'))
    x3min=float(abin.header('mesh','x3min'))
    x1max=float(abin.header('mesh','x1max'))
    x2max=float(abin.header('mesh','x2max'))
    x3max=float(abin.header('mesh','x3max'))
    xyz=np.array([x1min,x1max,x2min,x2max,x3min,x3max])
    # load data into yt
    dens_arr=abin.get_data('dens',xyz=list(xyz*1),level=0).T #to account for Athena's zyx
    temp_arr=abin.get_data('temp',xyz=list(xyz*1),level=0).T
    velx_arr=abin.get_data('velx',xyz=list(xyz*1),level=0).T
    vely_arr=abin.get_data('vely',xyz=list(xyz*1),level=0).T
    velz_arr=abin.get_data('velz',xyz=list(xyz*1),level=0).T
    pres_arr=abin.get_data('pres',xyz=list(xyz*1),level=0).T
    data = dict(density = (dens_arr*unit.density_cgs, "g*cm**-3"), 
                number_density = (dens_arr, "cm**-3"), 
                H_nuclei_density = (mu/muH*dens_arr, "cm**-3"), 
                temperature= (temp_arr*unit.temperature_cgs, "K"),
                pressure= (pres_arr*unit.pressure_cgs, "erg*cm**-3"),            
                velocity_x = (velx_arr*unit.velocity_cgs, "cm*s**-1"),
                velocity_y = (vely_arr*unit.velocity_cgs, "cm*s**-1"),
                velocity_z = (velz_arr*unit.velocity_cgs, "cm*s**-1"),
            )
    bbox = np.array([[x1min,x1max], [x2min,x2max], [x3min,x3max]])
    ds = yt.load_uniform_grid(data, dens_arr.shape, 
                              sim_time = abin.time,
                              length_unit="pc", periodicity=(False,False,False),
                              bbox=bbox, nprocs=256,default_species_fields='ionized')
    return ds

if __name__ == "__main__":
    import os
    myid = int(os.environ["SLURM_ARRAY_TASK_ID"])
    # initialization
    basedir = '/home/mg9443/scratch/share/data/snr/snr_net/'
    outbase = '/scratch/gpfs/nn0933/'
    baseid = 'Blast.hydro_w.'
    if len(sys.argv) == 2:
        path=os.path.join(basedir,f'{sys.argv[1]}/data/bin/') #os.path.join
        outdir = os.path.join(outbase,sys.argv[1])
        os.makedirs(outdir,exist_ok=True)
    else:
        print('Wrong argument')
        sys.exit()
    
    #check number of files
    files= sorted(glob.glob(os.path.join(path,baseid + "*")))
    fstart= int(files[0][-6:-4])
    fend= int(files[-1][-6:-4])
    # ilist = list(range(fstart, fend + 1))
    myilist = list(range(myid,fend,10))
    # print("available snapshots:",ilist)
    print("my snapshots:",myilist)
    #ilist = [0,10]

    # set up bin information for PDFs
    bin_info = setup_bin_field_info(Nbin=ds.domain_dimensions[0])

    # setup source model
    source_model = pyxsim.CIESourceModel("spex", 0.05, 11.0, 1000, 1.0, binscale="log")
    for i in myilist:
        abin = AthenaBinary(path=path,num=i,version='0.2')
        abin.load_binary(os.path.join(path,f"{baseid}{i:05d}.bin"))
        abin.config()
        # atbs.read([i]) # read binaries from these specific files, check number of files 
        # atbs.config()

        print(f"processing {i}")
        ds = ytload(abin)
        # add extra fields
        add_yt_fields(ds)

        # create xray field
        xray_fields = source_model.make_source_fields(ds, 0.5, 7.0) 
    
        # saving slice info
        slc = ds.slice('z',0)
        sfrb = slc.to_frb((64,'pc'), (ds.domain_dimensions[0],ds.domain_dimensions[0]))
        store= dict()
        for k in ['H_nuclei_density','temperature']:
            store[k]=sfrb[('gas',k)]
        save_dict_to_hdf5(store,os.path.join(outdir,f'slice_nT_{i:04d}.hdf5'))
        
        # create PDFs
        xflist = ['nH','nH' ,'velz','r' ,'R' , 'x']
        yflist = ['T' ,'pok','T'   ,'Lx','Lx', 'y']
        for xf,yf in zip(xflist, yflist):
            make_one_PDF(ds,xf,yf,bin_info, wflist=['Lx','vol','mass','etot','eloss'], 
                         inum=i,save=True,outdir=outdir)
            
        