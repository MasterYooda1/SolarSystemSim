import time
from SolarSystem import to_state_vec, save_interval, interval, step_count
import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
import astropy.units as u
import matplotlib.pyplot as plt
# Test to check whether the orbit remains stable after a certain time

body_names = ["Sun", "Mercury", "Venus", "Earth", "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Rogue"]
#delay = 3 # The Amount of days you want to test data in, requires the interval to be 86400 (1 day)
#if not delay % save_interval == 0:
#    raise ValueError("The Number of days you want to test falls between a save interval, please change the delay or the save interval")
    
chosen_body = "Earth"
days_difference = 5 # From the Day the Simulation was Started
new_testing_date = Time("2025-12-6 12:00:00.0", scale="tdb") # Time to Check for Accuracy

if not type(chosen_body) is str:
    raise TypeError("Only String Names, like 'Earth' are allowed")
if not chosen_body in body_names:
    raise ValueError("That Body is not in our Simulation")

pos, vel = get_body_barycentric_posvel(chosen_body.lower(), new_testing_date, ephemeris="jpl")
sun_pos, sun_vel = get_body_barycentric_posvel("sun", new_testing_date, ephemeris="jpl")

dataset_position, dataset_velocity = to_state_vec(pos, vel, new_testing_date)
dataset_position, dataset_velocity = np.array(dataset_position, dtype=float) , np.array(dataset_velocity, dtype=float)

# Load in the Test Data to Check for Accuracy
data_in = np.load("NBodyTestVerlet.npy", allow_pickle=True)

entry_no = (days_difference * 86400) // interval // save_interval # Value for the Entry we are testing. There are step_count entries, every interval seconds and saves every save_interval steps, converts days_difference to Seconds
if entry_no <= step_count:
    simulated_position = data_in[entry_no][4].position
    simulated_velocity = data_in[entry_no][4].velocity
else:
    raise IndexError("Out of Bounds! The data does not cover this far into the future, alter the interval and save_interval to fix this!")

change_in_position = simulated_position - dataset_position # Compute Inaccuracies
change_in_velocity = simulated_velocity - dataset_velocity

# Ugly One Liner to assign each variable with the relevant magnitudes 
dataset_pos_magnitude, simulated_pos_magnitude, dataset_vel_magnitude, simulated_vel_magnitude  = np.linalg.norm(dataset_position), np.linalg.norm(simulated_position), np.linalg.norm(dataset_velocity), np.linalg.norm(simulated_velocity)

pos_error = np.linalg.norm(change_in_position)
vel_error = np.linalg.norm(change_in_velocity)


# Checks to See if the values fall within a certain inaccuracy if not it throws an error
if abs(pos_error) > 10:
    raise ValueError(f"The Simulation is Not Accurate to {str(new_testing_date)}, the Position is more than 10km Off: {np.linalg.norm(change_in_position)/1000} km.")
if abs(vel_error) > 0.02:
    raise ValueError(f"The Simulation is Not Accurate to {str(new_testing_date)}, the Velocity is more than 0.02km/s Off: {np.linalg.norm(change_in_velocity)/1000} km/s.")








