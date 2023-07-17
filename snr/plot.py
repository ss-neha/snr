import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from athena_data_2023_01 import load_dict_from_hdf5, save_dict_to_hdf5

with open(f'./pdf_nH_pok_0000.p','rb') as fp:
    pdf_load = pickle.load(fp)
#pdf_loaded = pdf_nH_pok_0000.p
flist = ['vol','mass']
fig,axes = plt.subplots(1,5,figsize=(15,5))
for ax, wf in zip(axes,flist):
    plt.sca(ax)
    plt.pcolormesh(pdf_load['x_edges'],pdf_load['y_edges'],pdf_load[wf].T,
                   norm=LogNorm())
    plt.xlabel('Cylindrical Radius')
    plt.ylabel('XRay Luminosity')
    plt.yscale('log')
    #plt.xlim(right=20)
    plt.colorbar(label=wf)
plt.tight_layout()
plt.savefig(f'')

fig,axes = plt.subplots(1,5,figsize=(15,5))
#nH_edges = pdf_load['x_edges'][1]-pdf_load['x_edges'][0]
for ax, wf in zip(axes,flist):
    plt.sca(ax)
    plt.step(pdf_load['x_edges'][:-1],pdf[wf].sum(axis=1))
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(1.e-5,10)
    plt.xlabel('H Nuclei Density')
    plt.ylabel(wf)
    #plt.xlim(-1,20)
plt.tight_layout()

fig,axes = plt.subplots(1,5,figsize=(15,5))
#nH_edges = pdf_load['x_edges'][1]-pdf_load['x_edges'][0]
for ax, wf in zip(axes,flist):
    plt.sca(ax)
    plt.step(pdf_load['y_edges'][:-1],pdf[wf].sum(axis=0))
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(1.e-5,10)
    plt.xlabel('POK')
    plt.ylabel(wf)
    #plt.xlim(-1,20)
plt.tight_layout()

