# This file generates random initial conditions for the HVS     we simulate in the main file

import numpy as np
import random as rng
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


#
'''CONSTANTS AND PARAMETERS:'''
#

# amount of initial conditions generated
ini_con = 5000

# Plummer radius (core of the cluster) [kpc]
plum_r = 5

# transformation constant for km/s -> kpc/Myr
trafo = 1.023 / 1000

# baryonic mass of cluster [M_o]
mass_and = (0.5 / 6) * 4.211301617 * 10 ** 12

# maximum baryon density of cluster, calculated from Plummer model: rho(r = 0) [M_o / kpc^3]
rho_max = 6.702494692 * 10 ** 8

# maximum magnitude of velocity vector [kpc/Myr]
v_max = 1100 * trafo

# escape velocity of Andromeda as minimum magnitude of velocity vector [kpc/Myr]
v_min = 900 * trafo


#
'''ARRAYS:'''
#

# we generate the initial conditions in spherical coordinates at first

# define arrays for initial v, r, theta, phi and t values
# as well as v and r arrays for cartesian coordinates later

# array for position r-coordinate
r_coord = [None] * ini_con

# array for position azimuthal coordinate
phi_coord = [None] * len(r_coord)

# array for position inclination coordinate
theta_coord = [None] * len(r_coord)

# array for velocity magnitude
mag_vel = [None] * len(r_coord)

# array for velocity azimuthal coordinate
phi_vel = [None] * len(r_coord)

# array for velocity inclination coordinate
theta_vel = [None] * len(r_coord)

# array for initial time
t_coord = [None] * len(r_coord)

# arrays for Cartesian vector components
r1 = [None] * len(r_coord)

r2 = [None] * len(r_coord)

r3 = [None] * len(r_coord)

v1 = [None] * len(r_coord)

v2 = [None] * len(r_coord)

v3 = [None] * len(r_coord)


#
'''CALCULATIONS:'''
#

# define function for Plummer radius depending on density


def radius(a, rho, M):

    r = a * ((((4 * np.pi * a ** 3 * rho)/(3 * M)) ** (-0.4)) - 1) ** 0.5
    return r


# generate random floats in [0, 1] and convert them into
# random weighted radii, random theta and phi coordinates and random velocities

for i in range(len(r_coord)):

    # generate random numbers between 0-1
    temp1 = rng.random()
    # converting into random densities between 0 and rho_max
    dens = temp1 * rho_max

    # print(temp1)

    # convert them into random but weighted radii by inserting into the radius function
    x = radius(plum_r, dens, mass_and)
    r_coord[i] = x

    # isotropic angular coordinates
    # generate random velocity coordinates (spherical)
    phi_coord[i] = rng.uniform(0, 2 * np.pi)

    # random, but weighted theta coordinates / inclination angles
    theta_coord[i] = np.arccos(rng.uniform(-1, 1))

    # random send-off time between 10Gyrs and 12Gyrs
    t_coord[i] = rng.uniform(10000, 12000)


# separate loop for velocities
j = 0
# use while loop because we need to generate more magnitude than we use due to weighting
while j < len(r_coord):
    # random velocities between 900 and 1100, but high enough to escape Andromeda
    temp3 = rng.uniform(v_min, v_max)
    # weighting condition with function f(x)=exp(-12*x)
    temp4 = rng.uniform(0, np.exp(12 * v_min * (-1)))
    v = np.exp(12 * temp3 * (-1))
    # test if the generated random velocity fulfils criterium
    if v > temp4:
        # add velocity magnitude to array
        mag_vel[j] = temp3

        # generate isotropic angular coordinates for velocity
        phi_vel[j] = rng.uniform(0, 2 * np.pi)
        theta_vel[j] = np.arccos(rng.uniform(-1, 1))

        # for while loop
        j += 1

# list for test plot
# test_plot = [None] * len(r_coord)

# transform spherical coordinates to Cartesian
for i in range(len(r_coord)):
    r1[i] = r_coord[i] * np.sin(theta_coord[i]) * np.cos(phi_coord[i])

    r2[i] = r_coord[i] * np.sin(theta_coord[i]) * np.sin(phi_coord[i])

    r3[i] = r_coord[i] * np.cos(theta_coord[i])

    # repeat for velocities
    v1[i] = mag_vel[i] * np.sin(theta_vel[i]) * np.cos(phi_vel[i])

    v2[i] = mag_vel[i] * np.sin(theta_vel[i]) * np.sin(phi_vel[i])

    v3[i] = mag_vel[i] * np.cos(theta_vel[i])

    # test_plot[i] = mag_vel[i] / trafo


#
'''RESULTS:'''
#

# save all 6 initial conditions into a line in the textfile for main simulation file to read
np.savetxt("HVS-initial-conditions.txt", [])

# write column headers
f = open("HVS-initial-conditions.txt", "w")
f.write("initial_time" + " " + "x-position" + " " + "y-position" + " " + "z-position" + " " +
        "x-velocity" + " " + "y-velocity" + " " + "z-velocity" + " " + "velocity_magnitude" "\n")

# write in results
for i in range(len(r_coord)):
    # write in the cartesian coordinates for distance and velocity in each line
    f.write(str(t_coord[i]) + " " + str(r1[i]) + " " + str(r2[i]) + " " + str(r3[i]) + " " +
            str(v1[i]) + " " + str(v2[i]) + " " + str(v3[i]) + " " + str(mag_vel[i]) + "\n")


#
'''PLOTTING:'''
#

# optional plots for inital velocity, time and directions
'''
# plot distribution of initial velocities
plt.figure(1)
# plt.hist(test_plot, 30)
plt.hist(test_plot, bins=30,  weights=np.ones(len(test_plot)) / len(test_plot))
plt.xlabel('Initial velocity magnitude [km/s]')
plt.ylabel('% of test stars')
# plt.yscale('log')
# plt.xlim(0, 33)
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.savefig('initial_vel_mag.pdf')


# kartesian plot of star positions and velocities
X, Y, Z, U, V, W = zip(*test_plot[0:100])

plt.figure(2)
ax = plt.subplot(111, projection='3d')
ax.quiver(X, Y, Z, U, V, W)
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
ax.set_zlim([-5, 5])
ax.set_xlabel('x-axis [kpc]')
ax.set_ylabel('y-axis [kpc]')
ax.set_zlabel('z-axis [kpc]')
plt.show()


plt.hist(t_coord, 80)
plt.xlabel('Distance to center of Andromeda [kpc')
plt.ylabel('# of results')
plt.show()
'''
