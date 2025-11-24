import numpy as np
from matplotlib import pyplot as plt
from poliastro import constants
import astropy.units as u
# Find the Eccentricity of Any orbit
# To Get an accurate value you must scale the time interval with the moons period. (interval * 2000 = period_M)

GME = constants.GM_earth.to(u.m**3 / u.s**2).value

DataIn = np.load("NBodyTest.npy", allow_pickle=True)

# For Each time step
for entry in DataIn:
    # Exlude the Step Number from the List
    Bodies = entry[1:]
    # Find the Vector Difference between the Earth and Moon for Displacement and Velocity
    r = Bodies[4].position - Bodies[3].position
    v = Bodies[4].velocity - Bodies[3].velocity
    
    # Find the Scalar Magnitude of the distance between them
    rMag = np.linalg.norm(r)
    
    # Calculate the Angular Momentum
    h = np.cross(r, v)
    
    # Calculate the Vector Eccentricity of the Orbit
    # Uses the Equation e_vec = ((v x r) x r)/GM_earth - r^/|r| 
    e_vec = np.cross(v, h) / float(GME) - r/rMag
    
    # Calculate the Magnitude of the Eccentricity
    e = np.linalg.norm(e_vec)
    
# Find the Overall Mean Eccentricity
print(np.mean(e))
    
    
    