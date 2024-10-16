import numpy as np
import matplotlib.pyplot as plt

def BPTplots(exp,var = ['radius','velocity','meltRate'],amb = True):

    """
    plots outputs of BPTmodel
    From output dataset exp from BPTmodel, creates a figure with subplots of the variables specified in var. If var = 'all', will plot all 9 options.
    Depths of neutral density and maximum height will also be marked. If amb is True, ambient ocean profile will also be plotted with temperature, salinity, and density.

    See also BPTmodel

    Parameters
    ----------
    exp : xarray dataset
        result from a BPTmodel run 
    var : list of strings or 'all', optional
        list of strings containing variable names to plot. Variables may be any number of the following: radius, velocity, temperature, salinity, density, area, volumeFlux, momentumFlux, and meltRate.
        If set to 'all', will plot all 9 variable options.
    amb : boolean, optional
        select whether to include ambient ocean profile in appropriate plots

    Returns
    -------
    fig : matplotlib figure
    ax : list of axes output from plt.subplots(...)

    """

    fontsize = 16

    if var == "all":
        var = ['radius', 'velocity', 'temperature', 'salinity', 'density', 'area', 'volumeFlux', 'momentumFlux', 'meltRate']

    # determine number of rows and columns needed for number of var input
    if np.size(var)==9:
        m=3
        n=3
    else:
        m=int(np.floor(np.sqrt(np.size(var))))
        if np.remainder(np.size(var),m) == 0:
            n = int(np.size(var)/m)
        else:
            n=int(np.floor(np.size(var)/m)+1)


    fig, ax = plt.subplots(m,n,sharey=True,figsize=(10,10))

    for i in np.arange(np.size(var)):
        axt = ax.flatten()[i]
        axt.plot(exp['plume_'+var[i]],exp.depth,linewidth=2,label='Plume')
        if amb:
            if var[i] == 'temperature' or var[i] == 'salinity' or var[i] == 'density':
                axt.plot(exp['ambient_'+var[i]],exp.depth,linewidth=2,color='k',label='Ambient')
        
        axt.hlines(exp.plume_neutral_density_depth, 0, 1, transform=axt.get_yaxis_transform(),linestyle=':',label='Neutral Density')
        axt.hlines(exp.plume_maximum_height, 0, 1, transform=axt.get_yaxis_transform(),color='r',label='Maximum Height')
        if np.remainder(i,n)==0:
            axt.set_ylabel('depth (' + exp.depth.units + ')',fontsize=fontsize)
        else:
            axt.axes.get_yaxis().set_visible(False)
        axt.set_xlabel(var[i] + ' (' + exp['plume_'+var[i]].units + ')',fontsize=fontsize)
        if i == 0:
            axt.legend(loc = 'lower left',fontsize=fontsize)

    ax.flatten()[0].set_ylim(np.min(exp.depth),np.max(exp.depth))

    fig.tight_layout()

    return fig, ax



