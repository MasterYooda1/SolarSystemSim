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

def return_acc_at_pos(self, position, bodies):
    """Calculates the Acceleration at a position in space and returns it rather than updating the class like in compute_all_acceleration

    Args:
        position (np.array): position we are calculating from
        bodies (_type_): the bodies which we are calculating the acceleration based on

    Returns:
        _type_: _description_
    """
    acceleration = np.zeros(3)
    for body in bodies:
        if body is not self: 
            # Calculate the Vector Distance from one body to another
            dist_vec =  body.position - position
            # Calculate the Scalar Distance between the bodies
            dist = np.linalg.norm(dist_vec) 
            
            # Checks and Accounts for Division by Zero by replacing the 0 with a very small number
            min_dist = 1E-25
            if dist < min_dist:
                dist = min_dist
                
            # Unit Vector to Calculate the Acceleration along the right Axes 
            dir = dist_vec / dist
            
            # Adds the Acceleration from this body to the total acceleration of the body
            acceleration += np.array((self.G * body.mass/ dist ** 2) * dir, dtype=float)
    return acceleration

def update_verlet_velocity(bodies, interval):
    """Updates the Bodies using the Verlet Velocity algorithm, A Symplectic Integrator. Faster and more accurate. Energy Conserving.

    Args:
        bodies (list): List of particle objects containing celestial bodies
    """
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

def update_euler(bodies, interval):
    """Updates Position and Velocity using Euler's Method.

        Args:
            bodies (list): List of particle objects
            interval (Integer): Time Step
    """
    for body in bodies:
        body.compute_all_acceleration(body.position, bodies)
    
    for body in bodies:
        body.velocity += body.acceleration * interval
        body.position += body.velocity * interval
                
def update_RK4(bodies, interval):
    """A Non-Symplectic Integrator which uses multiple steps to refine the change in position and velocity, incredibly slow

    Args:
        bodies (list): list of particles
        interval (integer): time step
    """
    for body in bodies:
        pos_0 = body.position.copy()
        vel_0 = body.velocity.copy()
        
        # Step One, A Euler Approximation
        k1_vel = body.velocity.copy()
        k1_acc = return_acc_at_pos(body, pos_0, bodies)
        
        # Step Two, adds half the next step
        pos_1_mid = pos_0 + interval * k1_vel / 2
        k2_vel = vel_0 + interval * k1_acc / 2
        k2_acc = return_acc_at_pos(body, pos_1_mid, bodies)
        
        # Step Three
        pos_2_mid = pos_0 + interval * k2_vel / 2
        k3_vel = vel_0 + interval * k2_vel / 2
        k3_acc = return_acc_at_pos(body, pos_2_mid, bodies)
        
        # Step Four
        pos_3_mid = pos_0 + interval * k3_vel
        k4_vel = vel_0 + interval * k3_acc
        k4_acc = return_acc_at_pos(body, pos_3_mid, bodies)
        
        # Combine All the Steps and Update the Body
        body.position = pos_0 + (interval / 6) * (k1_vel + 2 * k2_vel + 2 * k3_vel + k4_vel)
        body.velocity = vel_0 + (interval / 6) * (k1_acc + 2 * k2_acc + 2 * k3_acc + k4_acc)

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
    

    # A Rogue Planet the Size of Earth, Orbitting around 1.3 times further than Pluto at 5x the speed
    #m_rogue = (constants.GM_earth / G).value
    #Rogue = Particle(
    #    position=np.array(position)*2,
    #    velocity=np.array(velocity) * np.array([-0.7, 1, 0]) * 1.5,
    #    acceleration=np.array([0, 0, 0]),
    #    name="Rogue",
    #    mass=m_rogue 
    #)
    #bodies.append(Rogue)
    
    return bodies

def run_sim(bodies, save_interval, step_count, interval):
    Data = [] # Data to be written to file

    for body in bodies:
        body.compute_all_acceleration(body.position, bodies)
    
    for step in range(step_count):
        # Loop through each body and update Positon and Velocity
        match chosen_updater:
            case 1:
                update_verlet_velocity(bodies, interval)
            case 2:
                update_euler(bodies, interval)
            case 3:
                update_RK4(bodies, interval)
        
        # If the Step is a Multiple of the save_interval Write to file.
        if step % save_interval == 0:
            Data.append([step, *(copy.deepcopy(b) for b in bodies)]) # Writes the Step number as well as a copy of the data for each body in the system


    np.save("NBodyTestVerlet.npy", Data, allow_pickle=True) # Writes the List Data to the aforementioned file
    
# Defining these variables outside of the main function allows them to be imported without running the code in main every time.
save_interval = 100 # Writes to File Every N Loops
interval = 50 # Updates the Simulation every N Seconds
step_count = 100000 # Repeats the update for N times
chosen_updater = 1 # The update method which the simulation uses (From 1 - 3)

def main():
    t = Time("2025-12-01 12:00:00", scale="tdb") # Desired Start Time,  1st December 2025
    bodies = loadBodies(t)

    run_sim(bodies, save_interval, step_count, interval)

if __name__ == "__main__": # Only Runs if You're not Exporting a Function to another file
    main()
    

    
    