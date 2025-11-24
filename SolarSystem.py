import numpy as np
from Particle import Particle
import copy
from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel
from poliastro import constants
from astropy.constants import G
from spiceypy import sxform, mxvg
import time

def to_state_vec(pos, vel, t):
    """
    Input:
    - Takes in the Cartesian Form of Astropy Data
    
    Function:
    - Converts it to usable np.array vector forms so we can do calculations with them
    
    Output:
    - The Position (np.array) and Velocity (np.array)
    """
    # Convert from default Units to the desired m and m/s and create a list of values as a vector.
    state_vec = [
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
    state_vececl = mxvg(trans, state_vec)
    
    # Assign the Converted Values to the needed vectors
    position = [state_vececl[0], state_vececl[1], state_vececl[2]]
    velocity = [state_vececl[3], state_vececl[4], state_vececl[5]]
    
    return position, velocity

def loadBodies(t):

    # Create the list of the Bodies in the Simulation
    bodies = []

    # The Sun
    # Get the Position and Velocity Vectors in Cartesian form from the JPL Ephemeris
    pos, vel = get_body_barycentric_posvel("sun", t, ephemeris="jpl")
    # Convert them to usable vectors
    position, velocity = to_state_vec(pos, vel, t)
    # Get the mass of the sun
    m_sun = (constants.GM_sun / G).value
    # Create an instance of Particle for the sun
    Sun = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Sun",
        mass=m_sun 
    )
    # Add the Sun to the list of stellar bodies
    bodies.append(Sun)

    # Mercury
    pos, vel = get_body_barycentric_posvel("mercury", t, ephemeris="jpl")
    position, velocity = to_state_vec(pos, vel, t)
    m_mercury = (constants.GM_mercury / G).value
    Mercury = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Mercury",
        mass=m_mercury 
    )
    bodies.append(Mercury)

    # Venus
    pos, vel = get_body_barycentric_posvel("venus", t, ephemeris="jpl")
    position, velocity = to_state_vec(pos, vel, t)
    m_venus = (constants.GM_venus / G).value
    Venus = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Venus",
        mass=m_venus 
    )
    bodies.append(Venus)

    # Earth
    pos, vel = get_body_barycentric_posvel("earth", t, ephemeris="jpl")
    position, velocity = to_state_vec(pos, vel, t)
    m_earth = (constants.GM_earth / G).value
    Earth = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Earth",
        mass=m_earth 
    )
    bodies.append(Earth)

    pos, vel = get_body_barycentric_posvel("moon", t, ephemeris="jpl")
    position, velocity = to_state_vec(pos, vel, t)
    m_moon = (constants.GM_moon / G).value
    Moon = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Moon",
        mass=m_moon 
    )
    bodies.append(Moon)

    # Mars
    pos, vel = get_body_barycentric_posvel("mars", t, ephemeris="jpl")
    position, velocity = to_state_vec(pos, vel, t)
    m_mars = (constants.GM_mars / G).value
    Mars = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Mars",
        mass=m_mars 
    )
    bodies.append(Mars)

    # Jupiter
    pos, vel = get_body_barycentric_posvel("jupiter", t, ephemeris="jpl")
    position, velocity = to_state_vec(pos, vel, t)
    m_jupiter = (constants.GM_jupiter / G).value
    Jupiter = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Jupiter",
        mass=m_jupiter 
    )
    bodies.append(Jupiter)

    # Saturn
    pos, vel = get_body_barycentric_posvel("saturn", t, ephemeris="jpl")
    position, velocity = to_state_vec(pos, vel, t)
    m_saturn = (constants.GM_saturn / G).value
    Saturn = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Saturn",
        mass=m_saturn 
    )
    bodies.append(Saturn)

    # Uranus
    pos, vel = get_body_barycentric_posvel("uranus", t, ephemeris="jpl")
    position, velocity = to_state_vec(pos, vel, t)
    m_uranus = (constants.GM_uranus / G).value
    Uranus = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Uranus",
        mass=m_uranus 
    )
    bodies.append(Uranus)

    # Neptune
    pos, vel = get_body_barycentric_posvel("neptune", t, ephemeris="jpl")
    position, velocity = to_state_vec(pos, vel, t)
    m_neptune = (constants.GM_neptune / G).value
    Neptune = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Neptune",
        mass=m_neptune 
    )
    bodies.append(Neptune)

    # Pluto
    pos, vel = get_body_barycentric_posvel("pluto", t, ephemeris="jpl")
    position, velocity = to_state_vec(pos, vel, t)
    m_pluto = (constants.GM_pluto / G).value
    Pluto = Particle(
        position=np.array(position),
        velocity=np.array(velocity),
        acceleration=np.array([0, 0, 0]),
        name="Pluto",
        mass=m_pluto 
    )
    bodies.append(Pluto)
    
    return bodies

def run_sim(bodies, save_interval, step_count, interval):
    Data = [] # Data to be written to file

    for body in bodies:
        body.compute_all_acceleration(bodies)
    
    for step in range(step_count):
        # Loop through each body and update Positon and Velocity
        for body in bodies:
            body.position += body.velocity * interval + 0.5 * body.acceleration * interval ** 2
            
        for body in bodies:
            body._next_acceleration = np.zeros(3, dtype=float)
            for body2 in bodies:
                if body is not body2:
                    dist_vec = body2.position - body.position
                    dist = np.linalg.norm(dist_vec)
                    if dist < 1E-25:
                        dist = 1E-25
                    direction = dist_vec / dist
                    body._next_acceleration += (body.G * body2.mass / dist ** 2) * direction
                    
        for body in bodies:
            body.velocity += 0.5 * (body.acceleration + body._next_acceleration) * interval
            body.acceleration = body._next_acceleration.copy()
        
        # If the Step is a Multiple of the save_interval Write to file.
        if step % save_interval == 0:
            Data.append([step, *(copy.deepcopy(b) for b in bodies)]) # Writes the Step number as well as a copy of the data for each body in the system


    np.save("NBodyTest.npy", Data, allow_pickle=True) # Writes the List Data to the aforementioned file
    
# Defining these variables outside of the main function allows them to be imported without running the code in main every time.
save_interval = 1 # Writes to File Every N Loops
interval = 10 # Updates the Simulation every N Seconds
step_count = 100000 # Repeats the update for N times

def main():
    t = Time("2025-11-17 14:25:00.0", scale="tdb") # Current Time
    bodies = loadBodies(t)

    run_sim(bodies, save_interval, step_count, interval)

if __name__ == "__main__": # Only Runs if You're not Exporting a Function to another file
    main()
    

    
    