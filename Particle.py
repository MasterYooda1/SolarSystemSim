import numpy as np
import math

class Particle():
    def __init__(
        self,
        position=np.array([0, 0, 0], dtype=float),
        velocity=np.array([0, 0, 0], dtype=float),
        acceleration=np.array([0, -10, 0], dtype=float),
        name='Ball',
        mass=1.0):
        
        
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.acceleration = np.array(acceleration, dtype=float)
        self.name = name
        self.mass = mass
        self.G = 6.67408E-11
        
    def __str__(self):
        return "Particle: {0}, Mass: {1:.3e}, Position: {2}, Velocity: {3}, Acceleration: {4}".format(
            self.name, self.mass,self.position, self.velocity, self.acceleration
        )
        
    def update(self, deltaT) -> None:
        self.position += self.velocity * deltaT
        self.velocity += self.acceleration * deltaT
        
    def getScalDistanceToBody(self, body):
        # Calculate the difference in the positions and use pythagorean theorem to find the scalar distance
        deltaPos = np.array(self.position - body.position, dtype=float)
        sum = 0
        for value in deltaPos:
            sum += value ** 2
            
        return math.sqrt(sum)
    
    def getVecDistanceToBody(self, body):
        # Calculate the difference in the positions
        return np.array(self.position - body.position, dtype=float)
            
    def updateGravitationalAcceleration(self, body):
        """
        Input:
        The Body in which the particle is interacting with (Assuming the body is a particle)
        
        Function:
        Update the acceleration of the particle based on the influence of the body.
        
        Output:
        None
        """
        # F_p = (G*M_p * M_b) / r^2
        # A_p = F_p/M_p = G * M_b / r^2
        dist_vec = self.getVecDistanceToBody(body)
        dist = self.getScalDistanceToBody(body)
        
        min_dist = 1E-25
        if dist < min_dist:
            dist = min_dist
            
        # Unit Vector
        dir = np.array(dist_vec / dist, dtype=float)
        
        self.acceleration = np.array(- (self.G * body.mass/ dist ** 2) * dir, dtype=float)
        
    def kineticEnergy(self):
        """
        Function
        ----------------
        Gets the Kinetic Energy of the Particle using K = 1/2 * m * v ** 2
        """
        v = np.linalg.norm(self.velocity)
        return 1/2 * self.mass * v ** 2
        
        
