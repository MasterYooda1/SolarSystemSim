from SolarSystem import toStateVec, saveInterval, bodies
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
import astropy.units as u
# Test to check whether the orbit remains stable after a certain time


bodyNames = []
delay = 3 # The Amount of days you want to test data in, requires the interval to be 86400 (1 day)
if not delay % saveInterval == 0:
    raise ValueError("The Number of days you want to test falls between a save interval, please change the delay or the save interval")

for body in bodies:
    bodyNames.append(body.name)
    
chosenBody = "Earth"
t = Time(f"2025-11-17 14:25:00.0", scale="tdb") # Current Time

if not type(chosenBody) is str:
    raise TypeError("Only String Names, like 'Earth' are allowed")
if not chosenBody in bodyNames:
    raise ValueError("That Body is not in our Simulation")

pos, vel = get_body_barycentric_posvel((chosenBody).lower(), t, ephemeris="jpl")

position, velocity = toStateVec(pos, vel)
position, velocity = np.array(position, dtype=float) * 1000, np.array(velocity, dtype=float) * 1000/86400

DataIn = np.load("TwoBodyTest.npy", allow_pickle=True)


