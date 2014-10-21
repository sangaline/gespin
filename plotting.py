import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
from functools import lru_cache
from gespin.core import nucleon_identities

@lru_cache(maxsize=2)
def build_sphere(theta_steps=6, phi_steps=12):
    theta = np.linspace(0, np.pi, num=theta_steps)
    phi = np.linspace(0, 2*np.pi, num=phi_steps)
    theta, phi = np.meshgrid(theta, phi)
    theta, phi = theta.ravel(), phi.ravel()
    mesh_x, mesh_y = (np.cos(phi)*theta, np.sin(phi)*theta)
    triangles = mpl.tri.Triangulation(mesh_x, mesh_y).triangles
    x, y, z = np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)
    return x, y, z, triangles

color_map = { nucleon_identities.unspecified : (208/255, 177/255, 33/255), \
              nucleon_identities.proton : (93/255, 46/255, 255/255), \
              nucleon_identities.neutron : (153/255, 90/255, 24/255), \
              nucleon_identities.antineutron : (196/255, 15/255, 96/255), \
              nucleon_identities.antiproton : (255/255, 56/255, 46/255) \
            }

def draw_nucleons(nucleons, alpha = 0.1):
    if not hasattr(nucleons, '__contains__'):
        return draw_nucleons([nucleons])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y, z, triangles = build_sphere()
    for nucleon in nucleons:
        #construct a sphere
        x_prime = nucleon.radius*x + nucleon.x
        y_prime = nucleon.radius*y + nucleon.y
        z_prime = nucleon.radius*z + nucleon.z
        r, g, b = color_map[nucleon.identity]
        color = (r, g, b, alpha)
        edgecolor = (r, g, b, 1)
        ax.plot_trisurf(x_prime, y_prime, triangles, z_prime, color=color, edgecolor=edgecolor, linewidth=0.1)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    return fig, ax
