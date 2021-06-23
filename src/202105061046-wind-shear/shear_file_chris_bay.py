import meshio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, sys, shutil
# File paths and directories
filename    = '../inputfiles/results/LES/high-ti/shear/array.mean0D_UAvg.vtk'
# Read the VTK file using meshio
print('Reading VTK file using meshio.')
mesh = meshio.read(filename=filename,
                    file_format='vtk')
# Get mesh points
points = mesh.points
mesh_type = list(mesh.cells_dict.keys())[0]
cells = list(mesh.cells_dict.values())[0]
# Convert points and cells to cell centers
cell_corners = points[cells]
cell_centers = np.mean(cell_corners, axis=1)
U_data = mesh.point_data['UAvg']
points = mesh.points
# Save as feather
df = pd.DataFrame({
    'x': points[:, 0],
    'y': points[:, 1],
    'z': points[:, 2],
    'u': U_data[:, 0],
    'v': U_data[:, 1],
    'w': U_data[:, 2] 
})
df_plane = df[df['x'] == np.min(df['x'])]
avg = []
for val in np.unique(df_plane['z']):
    tmp_avg = np.mean(df_plane[df_plane['z'] == val]['u'])
    print(tmp_avg)
    avg.append(tmp_avg)
z_vals = np.unique(df_plane['z'])
shear_fit = (z_vals / 90.0) ** 0.10 * 8.1

np.savetxt("uaveraged.txt", np.c_[z_vals,avg])
plt.figure()
plt.plot(avg, z_vals, 'k', label='SOWFA')
plt.plot(shear_fit, z_vals, 'r', label='fit')
plt.legend()
plt.show()