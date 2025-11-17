import numpy as np
from Particle import Particle
import copy

earthMass = 5.97237e24     # https://en.wikipedia.org/wiki/Earth
earthRadius = 63710 * 1e3  # https://en.wikipedia.org/wiki/Earth
Earth = Particle(
    position=np.array([0, 0, 0]),
    velocity=np.array([0, 0, 0]),
    acceleration=np.array([0, 0, 0]),
    name="Earth",
    mass=earthMass
)


satPosition = earthRadius + (35786 * 1e3)
satVelocity = np.sqrt(Earth.G * Earth.mass / satPosition)  # from centrifugal force = gravitational force
Satellite = Particle(
    position=np.array([satPosition, 0, 0]),
    velocity=np.array([0, satVelocity, 0]),
    acceleration=np.array([0, 0, 0]),
    name="Satellite",
    mass=100.
)

stepCount = 200000
interval = 6
Data = []

for i in range(0, stepCount * interval, interval):
    for particle1 in [Earth, Satellite]:
        for particle2 in [Earth, Satellite]:
            if particle1 != particle2:
                particle1.updateGravitationalAcceleration(particle2)
                particle1.update(interval)
    if (i / 6) % 100 == 0:
        # Append Data to File
        Data.append([i + 6, copy.deepcopy(Earth), copy.deepcopy(Satellite)])

np.save("TwoBodyTest.npy", Data, allow_pickle=True)
                
