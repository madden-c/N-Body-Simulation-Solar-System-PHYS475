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
Now we create the CelestialBody class. This is a fundamental component of the code that allows us to store the physical properties of the planets we are looking to simulate.
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
```
