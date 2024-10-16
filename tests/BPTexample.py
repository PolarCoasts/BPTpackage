
import numpy as np
from BPTmodel.model import BPTmodel
from BPTmodel.plots import BPTplots
import matplotlib.pyplot as plt
from BPTmodel.melt import melt_calc

## This example uses an arbitrary ambient profile (cold/fresh at surface, warm/salty at depth)
ambD=np.arange(-2,-300,-2)
ambT = 3-1/np.logspace(0.1,1.5,np.size(ambD))
ambS = 28-5/np.logspace(0.1,1.8,np.size(ambD))

## Define grounding line depth and discharge flux
gl=-290; # (m)
q=200;   # (m^3/s)

## Run model with default values and plot all plume properties
ex1=BPTmodel(ambD,ambT,ambS,gl,q)

fig1,ax1=BPTplots(ex1,'all')
fig1.suptitle('Line Plume')
fig1.tight_layout()

fig1.savefig('plots/LinePlumeExample.png')
plt.close(fig1)

## Run model as a point plume, specify entrainment coefficient, and plot just radius and vertical velocity
ex2=BPTmodel(ambD,ambT,ambS,gl,q,type='point',alpha=.08)

fig2,ax2=BPTplots(ex2,["radius", "velocity"])
fig2.suptitle('Point Plume')
fig2.tight_layout()

fig2.savefig('plots/PointPlumeExample.png')
plt.close(fig2)

## Run model for stacked ambient melt plumes and plot default properties
noq=1e-10; #small discharge to initiate plume
uhoriz=0.2; #ambient along-terminus velocity

ex3=BPTmodel(ambD,ambT,ambS,gl,noq,uh=uhoriz,type='stacked')

fig3,ax3=BPTplots(ex3)
fig3.suptitle('Stacked Ambient Melt Plumes')
fig3.tight_layout()

fig3.savefig('plots/StackedPlumeExample.png')
plt.close(fig3)

## Run independent melt calculation
# point estimate of melt
melt, tb, _ =melt_calc(-20,uhoriz,4,30)

# profile estimate of melt, specifying ice temp
uprof=np.linspace(.3,.01,np.size(ambD))
meltprof, _, _=melt_calc(ambD,uprof,ambT,ambS,iceT=0)
# convert to daily melt rate
meltprof=meltprof*24*60*60

# plot melt profile
figm, (axu,axm) = plt.subplots(1,2,sharey=True)
axu.plot(uprof,ambD,linewidth=2,color='k')
axu.set_ylabel('Depth (m)',fontsize=16)
axu.set_xlabel('Along-ice Velocity (m/s)',fontsize=16)
axm.plot(meltprof,ambD,linewidth=2)
axm.axes.get_yaxis().set_visible('False')
axm.set_xlabel('Melt Rate (m/day)',fontsize=16)
figm.suptitle('Ambient Melt - No Plume',fontsize=16)
figm.tight_layout()

figm.savefig('plots/MeltProfileExample.png')
plt.close(figm)


