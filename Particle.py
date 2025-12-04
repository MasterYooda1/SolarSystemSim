import numpy as np

"""
Potential integrators that can be Used:
IAS15 - https://arxiv.org/abs/1409.4779, https://academic.oup.com/mnras/article/446/2/1424/2892331
RK4
Symplectic Integrators
"""

class Particle():
    
    def __init__(
        self,
        position,
        velocity,
        acceleration,
        name,
        mass):
        """Constructor for Particle Class

        Args:
            position (np.array(float), optional): Position of the Particle in 3D space.
            velocity (np.array(float), optional): Velocity of the Particle in 3D Space.
            acceleration (np.array(float), optional): Acceleration of the Particle in 3D Space.
            name (str, optional): Name of the Particle.
            mass (float, optional): Mass of the Particle. Defaults to 1.0.
        """
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.acceleration = np.array(acceleration, dtype=float)
        self.name = name
        self.mass = mass
        self.G = 6.67408E-11
        self._next_acceleration = np.zeros(3, dtype=float)
        
    def __str__(self):
        """Returns the Values for the particle in string form

        Returns:
            String: a string containing the particles information
        """
        return "Particle: {0}, Mass: {1:.3e}, Position: {2}, Velocity: {3}, Acceleration: {4}".format(
            self.name, self.mass,self.position, self.velocity, self.acceleration
        )
        
    def update_euler(self, bodies, dt) -> None:
        """Updates the position and velocity of the particle. Very rudimentary, low intensity.

        Args:
            dt (Integer): the elapsed time since the last update
        """
        self.acceleration = self.compute_all_acceleration(self.position, bodies)
        self.velocity += self.acceleration * dt
        self.position += self.velocity * dt
    
    def update_verlet_velocity(self, dt):
        """Updates Position and Velocity using Verlet-Velocity Integration. Slightly more advanced and accurate than Euler, without being much more intense. !Be Cautious with Acceleration Updates

        Args:
            particles (list): List of Particle Objects
            dt (Integer): Time Step

        Returns:
            np.array(float): returns position and velocity vector(float)
        """
        self.position += self.velocity * dt + 0.5 * self.acceleration * dt ** 2
        
        self.velocity += 0.5 * (self.acceleration + self._next_acceleration) * dt
        
        self.acceleration = self._next_acceleration.copy()
        
    def update_RK4(self, particles, dt):
        #! NOT UPDATED TO NEW SYSTEM
        """Updates Position and Velocity using Runge-Kutta 4, the function denotes _ as a subscript not snake case.
        A lot more complicated, seems to either be bugged or the method itself is not as good, yields the same accuracy as Verlet-Velocity as per the Stability Test and Energy Test whilst taking twice as long as Verlet-Velocity.
        Args:
            particles (list): List of all the particles in the simulation
            dt (Integer): Time Step
        """
        pos_0 = self.position.copy() # Make a copy rather than setting equal to so its not changed as the position changes.
        vel_0 = self.velocity.copy()
        
        # First, Euler Approximation
        k1_vel = self.velocity.copy()
        k1_acc = self.compute_all_acceleration(pos_0, particles)
        
        pos_1_mid = pos_0 + dt * k1_vel / 2
        vel_1_mid = vel_0 + dt * k1_acc / 2
        k2_vel = vel_1_mid
        k2_acc = self.compute_all_acceleration(pos_1_mid, particles)
        
        pos_2_mid = pos_0 + dt * k2_vel / 2
        vel_2_mid = vel_0 + dt * k2_acc / 2
        k3_vel = vel_2_mid
        k3_acc = self.compute_all_acceleration(pos_2_mid, particles)
        
        pos_3_mid = pos_0 + dt * k3_vel
        vel_3_mid = vel_0 + dt * k3_acc
        k4_vel = vel_3_mid
        k4_acc = self.compute_all_acceleration(pos_3_mid, particles)
        
        pos_new = pos_0 + (dt/6.0) * (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel)
        vel_new = vel_0 + (dt/6.0) * (k1_acc + 2*k2_acc + 2*k3_acc + k4_acc)
        
        return pos_new, vel_new   
        
    def compute_all_acceleration(self, position, particles):
        acceleration = np.zeros(3)
        for body in particles:
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
        self.acceleration = acceleration
        
    def kinetic_energy(self):
        """Calculates the body's Kinetic Energy

        Returns:
            Float: Value of K in Joules as a Float
        """
        v = np.linalg.norm(self.velocity) # Calculates the Magnitude of the Velocity Vector
        return 1/2 * self.mass * v ** 2 # Returns the Kinetic Energy
        
        
