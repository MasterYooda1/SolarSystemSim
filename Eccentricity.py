import numpy as np
from matplotlib import pyplot as plt
from poliastro import constants
import astropy.units as u
# Find the Eccentricity of Any orbit
# To Get an accurate value you must scale the time interval with the moons period. (interval * 2000 = period_M)

gm_e = constants.gm_earth.to(u.m**3 / u.s**2).value

data_in = np.load("NBodyTest.npy", allow_pickle=True)

# For Each time step
for entry in data_in:
    # Exlude the Step Number from the List
    bodies = entry[1:]
    # Find the Vector Difference between the Earth and Moon for Displacement and Velocity
    r = bodies[4].position - bodies[3].position
    v = bodies[4].velocity - bodies[3].velocity
    
    # Find the Scalar Magnitude of the distance between them
    r_mag = np.linalg.norm(r)
    
    # Calculate the Angular Momentum
    h = np.cross(r, v)
    
    # Calculate the Vector Eccentricity of the Orbit
    # Uses the Equation e_vec = ((v x r) x r)/gm_earth - r^/|r| 
    e_vec = np.cross(v, h) / float(gm_e) - r/r_mag
    
    # Calculate the Magnitude of the Eccentricity
    e = np.linalg.norm(e_vec)
    
# Find the Overall Mean Eccentricity
print(np.mean(e))
    
    
    