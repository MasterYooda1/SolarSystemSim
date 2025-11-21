import numpy as np
from Particle import Particle
import copy
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
from poliastro import constants
from astropy.constants import G
from spiceypy import sxform, mxvg

def toStateVec(pos, vel, t):
    """
    Input:
    - Takes in the Cartesian Form of Astropy Data
    
    Function:
    - Converts it to usable np.array vector forms so we can do calculations with them
    
    Output:
    - The Position (np.array) and Velocity (np.array)
    """
    # Convert from default Units to the desired m and m/s and create a list of values as a vector.
    statevec = [
        pos.xyz[0].to("m").value,
        pos.xyz[1].to("m").value,
        pos.xyz[2].to("m").value,
        vel.xyz[0].to("m/s").value,
        vel.xyz[1].to("m/s").value,
        vel.xyz[2].to("m/s").value,
    ]
    
    # Create a tranform from an equitorial to ecliptic plane
    trans = sxform("J2000", "ECLIPJ2000", t.jd)
    
    # Perform the Transform using a Vector Matrix Multiplication
    statevececl = mxvg(trans, statevec)
    
    # Assign the Converted Values to the needed vectors
    position = [statevececl[0], statevececl[1], statevececl[2]]
    velocity = [statevececl[3], statevececl[4], statevececl[5]]
    
    return position, velocity

def loadBodies(t):

    # Create the list of the Bodies in the Simulation
    bodies = []

    # The Sun
    # Get the Position and Velocity Vectors in Cartesian form from the JPL Ephemeris
    pos, vel = get_body_barycentric_posvel("sun", t, ephemeris="jpl")
    # Convert them to usable vectors
    position, velocity = toStateVec(pos, vel, t)
    # Get the mass of the sun
    mSun = (constants.GM_sun / G).value
    # Create an instance of Particle for the sun
    Sun = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Sun",
        mass=mSun 
    )
    # Add the Sun to the list of stellar bodies
    bodies.append(Sun)

    # Mercury
    pos, vel = get_body_barycentric_posvel("mercury", t, ephemeris="jpl")
    position, velocity = toStateVec(pos, vel, t)
    mMercury = (constants.GM_mercury / G).value
    Mercury = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Mercury",
        mass=mMercury 
    )
    bodies.append(Mercury)

    # Venus
    pos, vel = get_body_barycentric_posvel("venus", t, ephemeris="jpl")
    position, velocity = toStateVec(pos, vel, t)
    mVenus = (constants.GM_venus / G).value
    Venus = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Venus",
        mass=mVenus 
    )
    bodies.append(Venus)

    # Earth
    pos, vel = get_body_barycentric_posvel("earth", t, ephemeris="jpl")
    position, velocity = toStateVec(pos, vel, t)
    mEarth = (constants.GM_earth / G).value
    Earth = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Earth",
        mass=mEarth 
    )
    bodies.append(Earth)

    pos, vel = get_body_barycentric_posvel("moon", t, ephemeris="jpl")
    position, velocity = toStateVec(pos, vel, t)
    mMoon = (constants.GM_moon / G).value
    Moon = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Moon",
        mass=mMoon 
    )
    bodies.append(Moon)

    # Mars
    pos, vel = get_body_barycentric_posvel("mars", t, ephemeris="jpl")
    position, velocity = toStateVec(pos, vel, t)
    mMars = (constants.GM_mars / G).value
    Mars = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Mars",
        mass=mMars 
    )
    bodies.append(Mars)

    # Jupiter
    pos, vel = get_body_barycentric_posvel("jupiter", t, ephemeris="jpl")
    position, velocity = toStateVec(pos, vel, t)
    mJupiter = (constants.GM_jupiter / G).value
    Jupiter = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Jupiter",
        mass=mJupiter 
    )
    bodies.append(Jupiter)

    # Saturn
    pos, vel = get_body_barycentric_posvel("saturn", t, ephemeris="jpl")
    position, velocity = toStateVec(pos, vel, t)
    mSaturn = (constants.GM_saturn / G).value
    Saturn = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Saturn",
        mass=mSaturn 
    )
    bodies.append(Saturn)

    # Uranus
    pos, vel = get_body_barycentric_posvel("uranus", t, ephemeris="jpl")
    position, velocity = toStateVec(pos, vel, t)
    mUranus = (constants.GM_uranus / G).value
    Uranus = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Uranus",
        mass=mUranus 
    )
    bodies.append(Uranus)

    # Neptune
    pos, vel = get_body_barycentric_posvel("neptune", t, ephemeris="jpl")
    position, velocity = toStateVec(pos, vel, t)
    mNeptune = (constants.GM_neptune / G).value
    Neptune = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Neptune",
        mass=mNeptune 
    )
    bodies.append(Neptune)

    # Pluto
    pos, vel = get_body_barycentric_posvel("pluto", t, ephemeris="jpl")
    position, velocity = toStateVec(pos, vel, t)
    mPluto = (constants.GM_pluto / G).value
    Pluto = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Pluto",
        mass=mPluto 
    )
    bodies.append(Pluto)
    
    return bodies

def runSim(bodies, saveInterval, stepCount, interval):
    Data = [] # Data to be written to file

    for step in range(stepCount):
        # Set the Magnitude of Acceleration to 0 for each body so it can be recalculated every step
        for body in bodies:
            body.acceleration = np.array([0.0, 0.0, 0.0], dtype=float)

        # Calculate the Acceleration for Body 2 on Body 1 and Sum to the overall acceleration vector
        # Has a big o of O(n^2) so for a large quantity of bodies it will be inefficient.
        for i, body1 in enumerate(bodies):
            for j, body2 in enumerate(bodies):
                if i != j: # If the Body is not itself 
                    body1.updateGravitationalAcceleration(body2)

        # Loop through each body and update Positon and Velocity
        for body in bodies:
            body.update(interval)

        # If the Step is a Multiple of the saveInterval Write to file.
        if step % saveInterval == 0:
            Data.append([step, *(copy.deepcopy(b) for b in bodies)]) # Writes the Step number as well as a copy of the data for each body in the system

    np.save("TwoBodyTest.npy", Data, allow_pickle=True) # Writes the List Data to the aforementioned file
    
saveInterval = 1 # Writes to File Every N Loops
interval = 50 # Updates the Simulation every N Seconds
stepCount = 20000 # Repeats the update for N times

def main():
    t = Time("2025-11-17 14:25:00.0", scale="tdb") # Current Time
    bodies = loadBodies(t)

    runSim(bodies, saveInterval, stepCount, interval)

if __name__ == "__main__": # Only Runs if You're not Exporting a Function to another file
    main()
    

    
    