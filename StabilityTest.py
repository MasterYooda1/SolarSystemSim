import time
from SolarSystem import toStateVec, saveInterval, interval, stepCount
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
import astropy.units as u
import matplotlib.pyplot as plt
# Test to check whether the orbit remains stable after a certain time

bodyNames = ["Sun", "Mercury", "Venus", "Earth", "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
#delay = 3 # The Amount of days you want to test data in, requires the interval to be 86400 (1 day)
#if not delay % saveInterval == 0:
#    raise ValueError("The Number of days you want to test falls between a save interval, please change the delay or the save interval")
    
chosenBody = "Earth"
daysDifference = 5
newTestingDate = Time("2025-11-22 14:25:00.0", scale="tdb") # Time to Check for Accuracy

if not type(chosenBody) is str:
    raise TypeError("Only String Names, like 'Earth' are allowed")
if not chosenBody in bodyNames:
    raise ValueError("That Body is not in our Simulation")

pos, vel = get_body_barycentric_posvel(chosenBody.lower(), newTestingDate, ephemeris="jpl")

datasetPosition, datasetVelocity = toStateVec(pos, vel, newTestingDate)
datasetPosition, datasetVelocity = np.array(datasetPosition, dtype=float) , np.array(datasetVelocity, dtype=float)

# Load in the Test Data to Check for Accuracy
DataIn = np.load("TwoBodyTest.npy", allow_pickle=True)

EntryNo = (daysDifference * 86400) // interval // saveInterval # Value for the Entry we are testing. There are stepCount entries, every interval seconds and saves every saveInterval steps, converts daysDifference to Seconds
print(EntryNo)
if EntryNo <= stepCount:
    simulatedPosition = DataIn[EntryNo][4].position
    simulatedVelocity = DataIn[EntryNo][4].velocity
else:
    raise IndexError("Out of Bounds! The data does not cover this far into the future, alter the interval and saveInterval to fix this!")

changeInPosition = DataIn[EntryNo][4].position - datasetPosition
changeInVelocity = DataIn[EntryNo][4].velocity - datasetVelocity

fig, ax = plt.subplots()

# Ugly One Liner
datasetPosMagnitude, simulatedPosMagnitude, datasetVelMagnitude, simulatedVelMagnitude  = np.linalg.norm(datasetPosition), np.linalg.norm(simulatedPosition), np.linalg.norm(datasetVelocity), np.linalg.norm(simulatedVelocity)

print(f"Dataset Position Magnitude: {datasetPosMagnitude} m")
print(f"Simulated Position Magnitude: {simulatedPosMagnitude} m")
print(f"Dataset Position Magnitude: {datasetVelMagnitude} m/s")
print(f"Dataset Position Magnitude: {simulatedVelMagnitude} m/s")

ax.plot(simulatedPosition[0], simulatedPosition[1], marker='o', label="Sim Position")
ax.plot(datasetPosition[0], datasetPosition[1], marker='o', label="JPL Position")
ax.plot(0, 0, marker='o', label='Sun')
ax.plot(1.496E11, 0, marker="o", label="Expected Distance")
ax.legend()


maxRange = 0
for vector in [datasetPosition, simulatedPosition]:
    maxRange = max(
        maxRange,
        np.max(np.abs(vector[0])),
        np.max(np.abs(vector[1])),
    ) + 1E11
    
ax.set_xlim(-maxRange, maxRange)
ax.set_ylim(-maxRange, maxRange)
    

plt.show()


if abs(np.linalg.norm(changeInPosition)) > 10000 or abs(np.linalg.norm(changeInVelocity)) > 10000:
    raise ValueError(f"The Simulation is Not Accurate to {str(newTestingDate)}, or the Position is more than 10km Off: {changeInPosition}, or The Simulation is Not Accurate to {str(newTestingDate)}, the Velocity is more than 10km/s Off: {changeInVelocity}")
else:
    print("Your values are in the expected range")








