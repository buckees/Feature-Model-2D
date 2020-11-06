"""Feature Model 2D. Mesh file."""
from FeatMod2d_ops import width, height, res_x, res_z
from FeatMod2d_mat import mat0, mat1, mat2, Si2d
from FeatMod2d_mesh import MESHGRID
import numpy as np
mesh = MESHGRID(width, height, res_x, res_z)
mesh.mat_input(Si2d)
mesh.mat = np.load('mesh_mat.npy')
mesh.find_surf()

for _idx, _surf in np.ndenumerate(mesh.surf):
    if _idx in mesh.surf_set:
        mesh.surf[_idx] = 2

mesh.plot()

# colMap = copy.copy(cm.get_cmap("Accent"))
# colMap.set_under(color='white')

# rec_surf = []
# for temp_idx in mesh.surf:
#     # temp_idx = (224, 13)
#     temp_svec, temp_stheta = mesh.calc_surf_norm(temp_idx, 
#                                                  radius=4,
#                                                  imode='Sum Vector')
#     rec_surf.append([temp_idx, temp_svec])

# def plot_surf_norm(ax, posn, svec):
#     ax.quiver(posn[0], posn[1],
#               svec[0], svec[1])


# fig, axes = plt.subplots(1, 2, figsize=(8, 8),
#                          constrained_layout=True)

# ax = axes[0]
# ax.contourf(mesh.x, mesh.z, mesh.mat, cmap=colMap, vmin=0.2, extend='both')
# ax.set_xlim(0.0, mesh.width)
# ax.set_ylim(0.0, mesh.height)

# ax = axes[1]
# ax.scatter(mesh.x, mesh.z, c=mesh.mat, s=1, cmap=colMap, vmin=0.2)
# ax.set_xlim(0.0, mesh.width)
# ax.set_ylim(0.0, mesh.height)

# for item in rec_surf:
#     temp_idx, temp_svec = item
#     temp_posn = np.array([mesh.x[temp_idx], mesh.z[temp_idx]])
#     plot_surf_norm(ax, temp_posn, temp_svec)

# plt.show()
# fig.savefig('init.png', dpi=600)
