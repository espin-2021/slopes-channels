import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from create_network import create_network_from_raster
from landlab.io.netcdf import (write_netcdf, read_netcdf)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from create_network import create_network_from_raster
from landlab.io.netcdf import (write_netcdf, read_netcdf)



filename = '/home/kpierce/Desktop/slopes-channels/space-exploring/presentation/simul/simul_00000.nc'
mg = read_netcdf(filename)


dx=100 # grid size
nmg = create_network_from_raster(
                                mg,
                                method='variable',
                                n_widths=dx,
                                min_channel_thresh=30000,
                                fields=['drainage_area', 'topographic__elevation']
                                )

# Shape of core nodes:
# to be modified for more complex grids (use a mask?)
number_of_core_node_rows = mg.number_of_node_rows - 2
number_of_core_node_columns = mg.number_of_node_columns - 2

#####################################################################

### Data for 3D plot of topographic surface
xplot = mg.node_x[mg.core_nodes].reshape((
        number_of_core_node_rows, number_of_core_node_columns))
yplot = mg.node_y[mg.core_nodes].reshape((
        number_of_core_node_rows, number_of_core_node_columns))

# Elevation of core nodes, reshaped for 3D plot:
zplot = mg.at_node['topographic__elevation'][mg.core_nodes].reshape(
        (number_of_core_node_rows, number_of_core_node_columns))

#####################################################################

### 3D plot of elevation surface:
# Figure and type of projection:
fig = plt.figure(1,tight_layout=True)
ax = plt.axes(projection='3d')

# Plot surface:
ax.plot_surface(xplot, yplot, zplot, cmap='pink', rstride=1, cstride=1,
                alpha=0.8)

# Set initial view of the graph (elevation and azimuth):
ax.view_init(elev=25, azim=-130)

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

ax.set_zlim(0,zplot.max()*2)

plt.axis('off')
#####################################################################

### also plot in stream traces...
#dvis=0.01*(zplot.max()-zplot.min())
#X,Y,Z = [],[],[]
#for link, nodes in enumerate(nmg.nodes_at_link):
#    x, y = nmg.x_of_node[nodes[0]], nmg.y_of_node[nodes[0]]
#    DX, DY = nmg.x_of_node[nodes[1]] - x, nmg.y_of_node[nodes[1]] - y
#    i, j = int(x/dx), int(y/dx)
#    z = zplot[i,j]
#    di = int(DX/dx)
#    dj=int(DY/dx)
#    DZ = zplot[i+di,j+dj]-z
#    ax.plot([x, x+DX], [y,y+DY],zs=[z+dvis,z+DZ+dvis],zorder=1,color='blue')
    

plt.show()
