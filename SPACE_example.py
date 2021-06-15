# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 13:35:03 2021

@author: laure
"""
#%%
"""
    Note: If timesteps are large enough that Es*dt (sediment erosion)
    exceeds sediment thickness H, the 'adaptive' solver is necessary to
    subdivide timesteps. Compare Es and H arrays to determine whether
    timesteps are appropriate or too large for the 'basic' solver.

"""
#%% IMPORT
import numpy as np
from landlab import RasterModelGrid, imshow_grid
from landlab.components import (FlowAccumulator,
                                DepressionFinderAndRouter,
                                Space,
                                FastscapeEroder)
import time

#%% Set up raster model grid
np.random.seed(seed = 5000)

mg = RasterModelGrid((20, 20), xy_spacing=100.0)
_ = mg.add_zeros('topographic__elevation', at='node')
mg.at_node['topographic__elevation'] += (mg.node_y / 10000. +
                                         mg.node_x / 10000. +
                                         np.random.rand(len(mg.node_y)) / 10.)
mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                       left_is_closed=True,
                                       right_is_closed=True,
                                       top_is_closed=True)
mg.set_watershed_boundary_condition_outlet_id(
    0, mg.at_node['topographic__elevation'], -9999.)
fsc_dt = 100.
space_dt = 100.

#%% Visualize grid

imshow_grid(mg, 'topographic__elevation')

#%% Instantiate Fastscape eroder, flow router, and depression finder

fr = FlowAccumulator(mg, flow_director='D8')
df = DepressionFinderAndRouter(mg)
fsc = FastscapeEroder(mg, K_sp=.001, m_sp=.5, n_sp=1)

#%% Burn in an initial drainage network using the Fastscape eroder:

for x in range(100):
    fr.run_one_step()
    df.map_depressions()
    fsc.run_one_step(dt=fsc_dt)
    mg.at_node['topographic__elevation'][0] -= 0.001 # Uplift

#%% Visualize grid

imshow_grid(mg, 'topographic__elevation')

#%% Add some soil to the drainage network:

_ = mg.add_zeros('soil__depth', at='node', dtype=float)
mg.at_node['soil__depth'] += 2.0
mg.at_node['topographic__elevation'] += mg.at_node['soil__depth']

#%% Visualize grid

imshow_grid(mg, 'topographic__elevation')

#%% Instantiate the Space component:

ha = Space(mg,
           K_sed=0.00001,
           K_br=0.00000000001,
           F_f=0.5,
           phi=0.1,
           H_star=1.,
           v_s=5.0,
           m_sp=0.5,
           n_sp = 1.0,
           sp_crit_sed=0,
           sp_crit_br=0,
           solver='basic')

#%% Now run the Space component for 2000 short timesteps:

start_time = time.time()

for x in range(2000):
    fr.run_one_step()
    df.map_depressions()
    ha.run_one_step(dt=space_dt)
    mg.at_node['bedrock__elevation'][0] -= 2e-6 * space_dt
    if x % 100 == 0:
        print("--", x*space_dt, "years -- %s seconds --" % round((time.time() - start_time), 1))

#%% Visualize grid

imshow_grid(mg, 'topographic__elevation')

#%%
#%%
#%%
"""
    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Shobe, C., Tucker, G., Barnhart, K. (2017). The SPACE 1.0 model: a Landlab
    component for 2-D calculation of sediment transport, bedrock erosion, and
    landscape evolution. Geoscientific Model Development  10(12), 4577 - 4604.
    https://dx.doi.org/10.5194/gmd-10-4577-2017

 """

#%%
"""Run step with CHILD-like solver that adjusts time steps to prevent slope
 flattening."""

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
import numpy as np

rg = RasterModelGrid((3, 4))
z = rg.add_zeros('topographic__elevation', at='node')
z[:] = 0.1 * rg.x_of_node
H = rg.add_zeros('soil__depth', at='node')
H += 0.1
br = rg.add_zeros('bedrock__elevation', at='node')
br[:] = z - H

fa = FlowAccumulator(rg, flow_director='FlowDirectorSteepest')
fa.run_one_step()
sp = Space(rg, K_sed=1.0, K_br=0.1,
           F_f=0.5, phi=0.0, H_star=1., v_s=1.0,
           m_sp=0.5, n_sp = 1.0, sp_crit_sed=0,
           sp_crit_br=0, solver='adaptive')
sp.run_one_step(dt=10.0)

np.round(sp.Es[5:7], 4)
array([ 0.0029,  0.0074])
np.round(sp.Er[5:7], 4)
array([ 0.0032,  0.0085])
np.round(H[5:7], 3)
array([ 0.088,  0.078])
