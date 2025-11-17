import numpy as np
import math

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
        
    def update(self, deltaT) -> None:
        """Updates the position and velocity of the particle

        Args:
            deltaT (Integer): the elapsed time since the last update
        """
        self.velocity += self.acceleration * deltaT
        self.position += self.velocity * deltaT
        
        
    def getScalDistanceToBody(self, body):
        """Calculates the Magnitude of the Vector between two bodies

        Args:
            body (Particle): The Body to find the Distance From

        Returns:
            float: the Scalar for the distance between two positions
        """
        # Calculate the difference in the positions and use pythagorean theorem to find the scalar distance
        deltaPos = np.array(self.position - body.position, dtype=float)
        sum = 0.0
        for value in deltaPos:
            sum += value ** 2
            
        return math.sqrt(sum)
    
    def getVecDistanceToBody(self, body):
        """Calculates the Difference in the 2 Position Vectors of the Bodies

        Args:
            body (Particle): The Body to find the Distance From

        Returns:
            np.array(float): the vector between the bodies
        """
        # Calculate the difference in the positions
        return np.array(body.position - self.position, dtype=float)
            
    def updateGravitationalAcceleration(self, body):
        """Computes the Gravitational Attraction between two bodies and adds it to the Current one

        Args:
            body (Particle): the Particle the current Particle is interacting with
        """
        
        # Calculate the Vector Distance from one body to another
        dist_vec = self.getVecDistanceToBody(body)
        # Calculate the Scalar Distance between the bodies
        dist = self.getScalDistanceToBody(body)
        
        # Checks and Accounts for Division by Zero by replacing the 0 with a very small number
        min_dist = 1E-25
        if dist < min_dist:
            dist = min_dist
            
        # Unit Vector to Calculate the Acceleration along the right Axes 
        dir = dist_vec / dist
        
        # Adds the Acceleration from this body to the total acceleration of the body
        self.acceleration += np.array((self.G * body.mass/ dist ** 2) * dir, dtype=float)
        
    def kineticEnergy(self):
        """Calculates the body's Kinetic Energy

        Returns:
            Float: Value of K in Joules as a Float
        """
        v = np.linalg.norm(self.velocity) # Calculates the Magnitude of the Velocity Vector
        return 1/2 * self.mass * v ** 2 # Returns the Kinetic Energy
        
        
