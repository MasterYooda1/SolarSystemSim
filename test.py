import numpy as np

DataIn = np.load("TwoBodyTest.npy", allow_pickle=True)

print("Testing reading the file TwoBodyTest.npy")

print("Printing Kinetic Energy of First Entry Earth and Satellite")
print("Earth KE: {:.5e}".format(DataIn[0][1].kineticEnergy()))
print("Satellite KE: {:.5e}".format(DataIn[0][2].kineticEnergy()))

print("Printing Kinetic Energy of Last Entry")
print("Earth KE: {:.5e}".format(DataIn[-1][1].kineticEnergy()))
print("Satellite KE: {:.5e}".format(DataIn[-1][2].kineticEnergy()))