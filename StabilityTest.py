import time
from SolarSystem import toStateVec, saveInterval
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
import astropy.units as u
# Test to check whether the orbit remains stable after a certain time

bodyNames = ["Sun", "Mercury", "Venus", "Earth", "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
#delay = 3 # The Amount of days you want to test data in, requires the interval to be 86400 (1 day)
#if not delay % saveInterval == 0:
#    raise ValueError("The Number of days you want to test falls between a save interval, please change the delay or the save interval")
    
chosenBody = "Earth"
newTestingDate = Time("2025-11-22 14:25:00.0", scale="tdb") # Time to Check for Accuracy

if not type(chosenBody) is str:
    raise TypeError("Only String Names, like 'Earth' are allowed")
if not chosenBody in bodyNames:
    raise ValueError("That Body is not in our Simulation")

pos, vel = get_body_barycentric_posvel(chosenBody.lower(), newTestingDate, ephemeris="jpl")

position, velocity = toStateVec(pos, vel, newTestingDate)
position, velocity = np.array(position, dtype=float) * 1000, np.array(velocity, dtype=float) * 1000/86400

# Load in the Test Data to Check for Accuracy
DataIn = np.load("TwoBodyTest.npy", allow_pickle=True)

changeInPosition = DataIn[8639][4].position * 1000 - position
changeInVelocity = DataIn[8639][4].velocity * 1000 - velocity

print(f"Difference in Position: {changeInPosition}")
print(f"Difference in Velocity: {changeInVelocity}")




