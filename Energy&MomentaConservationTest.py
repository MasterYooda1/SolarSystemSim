import numpy as np

# Load in the Test Data to Check for Accuracy
DataIn = np.load("TwoBodyTest.npy", allow_pickle=True)
bodies0 = DataIn[0][1:] # Take in Data from the first step

G = 6.67408E-11

# E_k = 1/2 * m * v^2
# U = -GMm/r
# Due to the Velocity the result would be a vector, since KE is a Scalar we must take the magnitude instead
E0_kinetic = 0
E0_potential = 0
for body in bodies0:
    E0_kinetic += 0.5 * body.mass * np.dot(body.velocity, body.velocity)

for i, body1 in enumerate(bodies0):
    for body2 in bodies0[i+1:]:
        r = np.linalg.norm(body2.position - body1.position)
        E0_potential -= G * body1.mass * body2.mass / r
        
E0_total = E0_kinetic + E0_potential

sampleStep = 1

bodiesFinal = DataIn[sampleStep][1:] # Take in the Data from the last step

# E_k = 1/2 * m * v^2
# U = -GMm/r
E_kinetic = 0
E_potential = 0
for body in bodiesFinal:
    E_kinetic += 0.5 * body.mass * np.dot(body.velocity, body.velocity)

for i, body1 in enumerate(bodiesFinal):
    for body2 in bodiesFinal[i+1:]:
        r = np.linalg.norm(body2.position - body1.position)
        E_potential -= G * body1.mass * body2.mass / r
        
energyTotal = E_kinetic + E_potential
print(f"Total Energy in the System: {energyTotal}")

energyLost = E0_total - energyTotal
print(f"Total Energy Lost in the System: {energyLost}")

energyPercentLost = (energyLost / energyTotal) * 100 
print(f"% Energy Lost: {abs(energyPercentLost)}%")
    
    
    
    
