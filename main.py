from Particle import Particle
import time as Time
import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as anim

def main():
    test1args = dict(
    position=numpy.array([0, 100, 0], dtype=float),
    velocity=numpy.array([0, 10, 0], dtype=float),
    acceleration=numpy.array([0, -10, 0], dtype=float),
    name="Ball",
    mass=500.0
    )
    Ball=Particle(**test1args)
    Ball.update(deltaT=0.1)
    animate_object(Ball)
        
def animate_object(object):
    fig = plt.figure(figsize=(5,3), dpi=200)
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel(r'x')
    ax.set_ylabel(r'y')
    ax.set_zlabel(r'z')
    ax.scatter(*object.position)
    
    plt.show()
    
    
    
if __name__ == "__main__":
    main()