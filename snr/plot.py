import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from athena_data_2023_01 import load_dict_from_hdf5, save_dict_to_hdf5

fig,axes = plt.subplots(1,5,figsize=(15,5))
for ax, wf in zip(axes,flist):
    plt.sca(ax)
    plt.pcolormesh(cr_edges,Lx_edges,pdf[wf].T,
                   norm=LogNorm())
    plt.xlabel('Cylindrical Radius')
    plt.ylabel('XRay Luminosity')
    plt.yscale('log')
    #plt.xlim(right=20)
    plt.colorbar(label=wf)
plt.tight_layout()

fig,axes = plt.subplots(1,5,figsize=(15,5))
for ax, wf in zip(axes,flist):
    plt.sca(ax)
    dR=np.diff(cr_edges)
    R = 0.5*(cr_edges[1:]+cr_edges[:-1])
    A = 2*np.pi*R*dR
    plt.step(cr_edges[:-1],pdf[wf].sum(axis=1)/A)
    #plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(1.e-5,10)
    plt.xlabel('Cylindrical Radius')
    plt.ylabel(wf)
    plt.xlim(-1,20)

plt.tight_layout()