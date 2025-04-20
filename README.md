# N-Body-Simulation-Solar-System-PHYS475
  In this repository, you will find the Python code for the N-body simulation of our solar sytstem's inner planetary bodies. The goal of this code is firstly, to visually simulate the gravitational interactions between inner solar system bodies in VPython, and secondly, to compare Euler's method and 4th order Runge-Kutta (RK4) integration techniques when used for N-body simulations. 

  We first start by importing our necessary libraries for this simulation. We will be using VPython, an open-source Python library that excels in imple and intuitive visualizations in three-dimensions. We will also be using AstroPy, a core package with numerous astronomy and astrophysics tools that will supply position values for us to measure our Euler and RK4 approximations against. We will also be using several other commonly used libraries such as pandas, numpy, matplotlib, and datetime.
```python
from vpython import *          # For 3D visualization
import pandas as pd            # For data handling and analysis
import numpy as np             # For numerical operations
from datetime import datetime, timedelta  # For time handling
import matplotlib.pyplot as plt  # For plotting results
from astropy.time import Time  # For astronomical time calculations
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric  # For planetary ephemeris
import astropy.units as u      # For unit handling in astronomy
```
Next, we define our physical constants and our simulation parameters. We define the gravitational constant (G), and we also include a radius_scale_factor to be optionally used to increase planetary bodies' visualized radius without altering the radii used for computations. We define our time step (dt), which will determine how often our simulation makes its calculations for position and velocity. The smaller we make dt, the more accurate our approximations should be. However, if we want to visualize multiple complete orbits, then we should increase dt, or total_simulation_time could be increased if you have the time to run a longer simulation with smaller time steps. 
```python
G = 6.674e-11                  # Gravitational constant (m^3 kg^-1 s^-2)
radius_scale_factor = 1e3      # Scaling factor for visualization
dt = 1e4                       # Time step (seconds)
total_simulation_time = dt * 1e4  # Total simulation duration 
bodies_to_compare = ['sun', 'mercury', 'venus', 'earth', 'mars']  # Bodies to simulate
```
Now we create the CelestialBody class. This is a fundamental component of the code that allows us to store the physical properties of the planets needed for the simulation. Once we establish the CelestialBody class, we have defined the structure of each body we look to simulate, now we provide the initial data for each object that will be used to perform the N-Body simulation. We define what visualization characteristics each body will have, and then use VPython vectors to hold initial position and velocity of each planetary body. Our data values for initial position and velocity are pulled from AstroPy JPL Ephemeris via the get_body_barycentric_posvel prompt. The position and velocity that are retrieved are in cartesian representation and are in AU and AU/day respectively. The position and velocity vectors which are hardcoded here were converted to the units of m and m/s respectively.
```python
class CelestialBody:
    """Class to represent celestial bodies with physical properties"""
    def __init__(self, name, color, pos, radius, mass, v):
        self.name = name       # Body name (e.g., "Earth")
        self.color = color     # Visualization color
        self.pos = pos         # Position vector (x,y,z) in meters
        self.radius = radius   # Physical radius in meters
        self.mass = mass       # Mass in kg
        self.v = v             # Velocity vector (vx,vy,vz) in m/s
        self.F = vector(0, 0, 0)  # Force vector (initially zero)
        self.p = self.v * self.mass  # Momentum vector (p = mv)

initial_bodies = [
    # Sun (central body)
    CelestialBody("Sun", color.yellow, vector(0,0,0), 
                 6.96342e8, 1.9884e30, vector(0,0,0)),
    
    # Mercury
    CelestialBody("Mercury", color.white, vector(-5.89406e10, -1.453e10, -6.1304e9), 
                 2.4397e6, 3.3011e23, vector(6591.68, -44823.7, -2404.34)),
    
    # Venus
    CelestialBody("Venus", color.orange, vector(9.33511e10, 5.48211e10, 4.42522e9), 
                 6.051e6, 4.5675e24, vector(-16959.9, 30567.2, 882.861)),
    
    # Earth
    CelestialBody("Earth", color.blue, vector(-1.38235e11, 6.09568e10, 0), 
                 6.378e6, 5.9735e24, vector(-12262.7, -26812.9, 0)),
    
    # Mars
    CelestialBody("Mars", color.red, vector(2.0567e11, -2.61614e10, 3.64974e9), 
                 1.794e6, 6.4171e24, vector(3621.45, 26168.2, 452.759)),
    
    # Jupiter
    CelestialBody("Jupiter", color.orange, vector(5.876e11, -3.885e11, -1.541e10),
                 6.9911e7, 1.89813e27, vector(6831.3, 10210.0, -191.5)),
    
    # Saturn
    CelestialBody("Saturn", color.yellow, vector(9.048e11, -1.075e12, -2.292e10),
                 5.8232e7, 5.6834e26, vector(6475.1, 5389.5, -397.3)),
    
    # Uranus
    CelestialBody("Uranus", color.cyan, vector(2.191e12, 1.964e12, -2.499e10),
                 2.5362e7, 8.6813e25, vector(-4574.3, 4991.7, 71.3)),
    
    # Neptune
    CelestialBody("Neptune", color.blue, vector(4.381e12, -8.28e11, -8.704e10),
                 2.4622e7, 1.0241e26, vector(966.8, 5495.1, -131.6))
]

```
