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
import matplotlib.pyplot as plt

from landlab import RasterModelGrid, imshow_grid
from landlab.components import (FlowAccumulator,
                                DepressionFinderAndRouter,
                                Space,
                                LinearDiffuser)
import time

#%% Set up parameters for model (units are m and yr)

dt = 50
runtime = 100000

uplift_rate = 0.0001
K_dif = 0.005

soil_depth_initial = 0.0
K_sed=0.0005
K_br=0.0001
F_f=0
phi=0
H_star=1.0
v_s=2.0
m_sp=0.5
n_sp = 1.0
sp_crit_sed=0
sp_crit_br=0
solver='basic'

#%% Set up grid size and spacing
num_rows = 30
num_cols = 30
dx = 100.

# #%% Set up time step imposed by CFL and von Neumann criteria
# dt_vN = dx*dx/(4.01*K_dif)  # von Neumann criterion for diffusion equation
# C_max = 1
# area = (num_cols*dx)*(num_rows*dx)
# dt_CFL = C_max*(dx/(K_sp*(area**m_sp)))   # CFL criterion for strea power equation
# dt = round(min(dt_vN,dt_CFL))

#%% Set up raster model grid
np.random.seed(seed = 5000)

mg = RasterModelGrid((num_rows, num_cols), xy_spacing=dx)
_ = mg.add_zeros('topographic__elevation', at='node')

# Random mm-scale noise
mg.at_node['topographic__elevation'] += (np.random.rand(len(mg.node_y)) / 100.)

# Sloped surface
mg.at_node['topographic__elevation'] += (mg.node_y / 100000. +
                                         mg.node_x / 100000.)

# Boundary conditions
mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                       left_is_closed=True,
                                       right_is_closed=True,
                                       top_is_closed=True)
mg.set_watershed_boundary_condition_outlet_id(
    0, mg.at_node['topographic__elevation'], -9999.)

# Soil depth
_ = mg.add_zeros('soil__depth', at='node', dtype=float)
mg.at_node['soil__depth'] += soil_depth_initial
mg.at_node['topographic__elevation'] += mg.at_node['soil__depth']

#%% Visualize grid
fig = plt.figure()
plot = plt.subplot()
imshow_grid(mg, 'topographic__elevation')
plt.show()

#%% Instantiate components:

fr = FlowAccumulator(mg, flow_director='D8')
df = DepressionFinderAndRouter(mg)
ld = LinearDiffuser(mg, linear_diffusivity=K_dif)
sp = Space(mg,
           K_sed=K_sed,
           K_br=K_br,
           F_f=F_f,
           phi=phi,
           H_star=H_star,
           v_s=v_s,
           m_sp=m_sp,
           n_sp=n_sp,
           sp_crit_sed=sp_crit_sed,
           sp_crit_br=sp_crit_br,
           solver=solver)

#%% Run model time loop:

start_time = time.time()

for x in range(runtime//dt):
    fr.run_one_step()
    df.map_depressions()
    sp.run_one_step(dt=dt)
    ld.run_one_step(dt=dt)
    mg.at_node['bedrock__elevation'][mg.core_nodes] += uplift_rate * dt
    if x*dt % 1000 == 0:
        print("--", x*dt, "years -- %s seconds --" % round((time.time() - start_time), 1))

#%% PLOT: topo map and sediment flux map

#Show topo map
fig = plt.figure()
plot = plt.subplot()
imshow_grid(mg, 'topographic__elevation', plot_name='Topographic Elevation', var_name = 'Topographic Elevation', var_units='m', grid_units=('m', 'm'), cmap='terrain')
plt.show()

#Show sediment flux map
fig = plt.figure()
plot = plt.subplot()
imshow_grid(mg, 'sediment__flux', plot_name='Sediment flux', var_name = 'Sediment flux', var_units=r'm$^3$/yr', grid_units=('m', 'm'), cmap='Blues')
plt.show()

#%% PLOT: slope-area

area = mg.at_node['drainage_area'][mg.core_nodes]
slope = mg.at_node['topographic__steepest_slope'][mg.core_nodes]

fig = plt.figure()
slope_area_plot = plt.subplot()
slope_area_plot.set_xscale('log')
slope_area_plot.set_yscale('log')
plt.scatter(area, slope)
plt.show()

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
