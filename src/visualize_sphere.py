"""Script to visualize the sphere object created for each atom"""

import matplotlib.pyplot as plt
from protein_structure import Sphere

# First we instanciate the sphere object 
atom_sphere = Sphere(radius=3.2, nbpoints=400)
atom_sphere.compute_points_coordinate()

# We extract the coordianetes from the sphere points
x = [pt[0] for pt in atom_sphere.pointlist]
y = [pt[1] for pt in atom_sphere.pointlist]
z = [pt[2] for pt in atom_sphere.pointlist]

# We create the plot
fig = plt.figure()  
ax = plt.axes(projection='3d')

# We plot the points
ax.scatter(x, y, z, color='blue')
ax.set_title('Point test on atom sphere surface')


ax.set_box_aspect([1, 1, 1])  

plt.show()