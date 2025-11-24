"""
Constants relevant to the simulation geometry and meshing.

"""

import numpy as np

# Parallelisation
num_cpus = 6
num_domains = num_cpus

# For sketching. Not important but needs specification
sheet_size = 10


# Rockwell indenter geometry specification
tip_radius = 0.2  # [mm]
cone_angle = 60  # [degrees]

# Geometry coordinates for Rockwell indenter -- Do not change!
xc1 = 0.0  # x coordinate of center of spherical tip
yc1 = 0.0  # y coordinate of center of spherical tip
xc2 = tip_radius * np.cos(-cone_angle * np.pi / 180)  # x coordinate of sphere end point
yc2 = tip_radius + tip_radius * np.sin(
    -cone_angle * np.pi / 180
)  # y coordinate of sphere end point
xc3 = tip_radius * np.cos(
    (-cone_angle - (90 - cone_angle) / 2.0) * np.pi / 180
)  # x coordinate of shpere mid point
yc3 = tip_radius + tip_radius * np.sin(
    (-cone_angle - (90 - cone_angle) / 2.0) * np.pi / 180
)  # y coordinate of sphere mid point
xl1 = xc2  # line x coorcinate 1
yl1 = yc2  # line y coordinate 1
xl2 = xl1 + 0.5 * np.cos((90 - cone_angle) * np.pi / 180)  # line x coordinate 2
yl2 = yl1 + 0.5 * np.sin((90 - cone_angle) * np.pi / 180)  # line y coordinate 2


# Substrate coordinates for 3D rectangle.
xs1 = 0.0  # x coordinate of first point
ys1 = 0.0  # y coordinate of first point
zs1 = 0.0  # z coordinate of first point
xs2 = 0.6  # x coordinate of second point
ys2 = 0.5  # y coordinate of second point
zs2 = 3.00  # z coordinate of second point - extrude depth


# For partitioning of substrate to allow better meshing
dpo_x = 0.24
dpo_y = 0.06
dpo_z = 0.25


# Scratch length and depth
scratch_length = 2
scratch_depth = -20e-3


# Coarse mesh sizes. Used on edge far from area of interest
coarse_mesh_size_0 = 0.05
coarse_mesh_size_1 = 0.15
coarse_mesh_size_2 = 0.30


# Analysis times
scratch_time = 0.01  # [s]
unload_time = 0.0001  # [s]


# Sample output frequencies
sample_frequency_unloading = unload_time / 20.0
sample_frequency_scratching = scratch_time / 20.0
sample_force_frequency = scratch_time / 100.0


# Max degredation of element stiffness
max_degradation = 0.6


# Naming of various Abaqus objects
indenter_name = "RockwellIndenter"
indenter_set_name = indenter_name + "Set"
indenter_instance_name = indenter_name + "Inst"

substrate_name = "Substrate"
substrate_set_name = substrate_name + "Set"
substrate_instance_name = substrate_name + "Inst"

master_surface_name = "m_Surf-1"
slave_surface_name = "s_Surf-1"
contact_region_nodes_name = "contactRegionNodes"
