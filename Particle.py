import numpy as np
import math

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
            position (np.array(float), optional): Position of the Particle in 3D space. Defaults to np.array([0, 0, 0], dtype=float).
            velocity (np.array(float), optional): Velocity of the Particle in 3D Space. Defaults to np.array([0, 0, 0], dtype=float).
            acceleration (np.array(float), optional): Acceleration of the Particle in 3D Space. Defaults to np.array([0, -10, 0], dtype=float).
            name (str, optional): Name of the Particle. Defaults to 'Ball'.
            mass (float, optional): Mass of the Particle. Defaults to 1.0.
        """
        
        
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.acceleration = np.array(acceleration, dtype=float)
        self.name = name
        self.mass = mass
        self.G = 6.67408E-11
        
    def __str__(self):
        """Returns the Values for the particle in string form

        Returns:
            String: a string containing the particles information
        """
        return "Particle: {0}, Mass: {1:.3e}, Position: {2}, Velocity: {3}, Acceleration: {4}".format(
            self.name, self.mass,self.position, self.velocity, self.acceleration
        )
        
    def updateEuler(self, bodies, dt) -> None:
        """Updates the position and velocity of the particle. Very rudimentary, low intensity.

        Args:
            dt (Integer): the elapsed time since the last update
        """
        self.acceleration = self.computeAllAcceleration(self.position, bodies)
        self.velocity += self.acceleration * dt
        self.position += self.velocity * dt
    
    def updateLeapfrog(self, particles, dt):
        """Updates Position and Velocity using Leapfrog Integration. Slightly more advanced and accurate than Euler, without being much more intense.

        Args:
            particles (list): List of Particle Objects
            dt (Integer): Time Step

        Returns:
            np.array(float): returns position and velocity vector(float)
        """
        pos0 = self.position.copy()
        vel0 = self.velocity.copy()
        
        aNow = self.computeAllAcceleration(pos0, particles)
        posNext = pos0 + vel0 * dt + 0.5 * aNow * dt ** 2
        
        aNext = self.computeAllAcceleration(posNext, particles)
        velNext = vel0 + 0.5 * (aNow + aNext) * dt
        
        return posNext, velNext
        
    
    def updateRK4(self, particles, dt):
        """Updates Position and Velocity using Runge-Kutta 4, the function denotes _ as a subscript not snake case.
        A lot more complicated, seems to either be bugged or the method itself is not as good, yields the same accuracy as Leapfrog as per the Stability Test and Energy Test whilst taking twice as long as Leapfrog.
        Args:
            particles (list): List of all the particles in the simulation
            dt (Integer): Time Step
        """
        pos0 = self.position.copy() # Make a copy rather than setting equal to so its not changed as the position changes.
        vel0 = self.velocity.copy()
        
        # First, Euler Approximation
        k1Vel = self.velocity.copy()
        k1Acc = self.computeAllAcceleration(pos0, particles)
        
        pos1Mid = pos0 + dt * k1Vel / 2
        vel1Mid = vel0 + dt * k1Acc / 2
        k2Vel = vel1Mid
        k2Acc = self.computeAllAcceleration(pos1Mid, particles)
        
        pos2Mid = pos0 + dt * k2Vel / 2
        vel2Mid = vel0 + dt * k2Acc / 2
        k3Vel = vel2Mid
        k3Acc = self.computeAllAcceleration(pos2Mid, particles)
        
        pos3Mid = pos0 + dt * k3Vel
        vel3Mid = vel0 + dt * k3Acc
        k4Vel = vel3Mid
        k4Acc = self.computeAllAcceleration(pos3Mid, particles)
        
        posNew = pos0 + (dt/6.0) * (k1Vel + 2*k2Vel + 2*k3Vel + k4Vel)
        velNew = vel0 + (dt/6.0) * (k1Acc + 2*k2Acc + 2*k3Acc + k4Acc)
        
        return posNew, velNew   
        
    def computeAllAcceleration(self, posMid, particles):
        acceleration = np.zeros(3)
        for body in particles:
            if body is not self: 
                # Calculate the Vector Distance from one body to another
                dist_vec =  body.position - posMid
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
        return acceleration.copy()
    
    def updateGravitationalAcceleration(self, body):
        """Computes the Gravitational Attraction between two bodies and adds it to the Current one

        Args:
            body (Particle): the Particle the current Particle is interacting with
        """
        
        
        
    def kineticEnergy(self):
        """Calculates the body's Kinetic Energy

        Returns:
            Float: Value of K in Joules as a Float
        """
        v = np.linalg.norm(self.velocity) # Calculates the Magnitude of the Velocity Vector
        return 1/2 * self.mass * v ** 2 # Returns the Kinetic Energy
        
        
