import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
M_earth = 5.972e24  # Earth's mass (kg)
M_moon = 7.342e22  # Moon's mass (kg)
R_earth = 6.371e6  # Earth's radius (m)
R_moon = 1.737e6  # Moon's radius (m)

# Function to calculate the gravitational force on a satellite including the Moon's effect
def gravitational_force_including_moon(mass_satellite, r):
    # Calculate the gravitational forces from both the Earth and the Moon
    force_earth = (G * M_earth * mass_satellite) / (np.linalg.norm(r) ** 2)
    distance_to_moon = np.linalg.norm(r - np.array([0, 0, 384400e3]))  # Distance to the Moon
    force_moon = (G * M_moon * mass_satellite) / (distance_to_moon ** 2)
    
    # Calculate the direction of the gravitational force from the Moon
    direction_to_moon = np.linalg.norm((np.array([0, 0, 384400e3]) - r) / distance_to_moon)
    
    # Calculate the net gravitational force including both the Earth and the Moon
    gravitational_force_total = force_earth - force_moon * direction_to_moon
    
    return gravitational_force_total

# Function to convert Cartesian coordinates to longitude, latitude, and altitude
def cartesian_to_geodetic(x, y, z):
    # WGS84 parameters
    a = 6378137.0  # Semi-major axis (m)
    b = 6356752.3142  # Semi-minor axis (m)
    
    # Calculate other params
    p = np.sqrt(x ** 2 + y ** 2)
    e2 = 1 - (b ** 2) / (a ** 2)
    
    ep = np.sqrt((a ** 2 - b ** 2) / (b ** 2))
    th = np.arctan2(a * z, b * p)
    
    # Calculate longitude
    longitude = np.arctan2(y, x)

    # Calculate latitude
    latA = z + ep ** 2 * b * np.sin(th) ** 3
    latB = p - e2 * a * np.cos(th) ** 3
    latitude = np.arctan2(latA, latB)
    
    # Calculate altitude
    N = a / np.sqrt(1 - e2 * np.sin(latitude) ** 2)
    altitude = p / np.cos(latitude) - N

    # Convert latitude and longitude from radians to degrees
    latitude = np.degrees(latitude)
    longitude = np.degrees(longitude)

    return longitude, latitude, altitude

# Function to calculate the satellite's position in orbit
def satellite_position(semi_major_axis, eccentricity, inclination, argument_of_perigee, true_anomaly):
    r = semi_major_axis * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos(true_anomaly))
    x = r * (np.cos(argument_of_perigee) * np.cos(true_anomaly) - np.sin(argument_of_perigee) * np.sin(true_anomaly) * np.cos(inclination))
    y = r * (np.sin(argument_of_perigee) * np.cos(true_anomaly) + np.cos(argument_of_perigee) * np.sin(true_anomaly) * np.cos(inclination))
    z = r * (np.sin(true_anomaly) * np.sin(inclination))
    return x, y, z

# Function to calculate satellite trajectories over time including the Moon's gravitational effect
def calculate_satellite_trajectories_including_moon(num_satellites, semi_major_axis, eccentricity, inclination, argument_of_perigee, masses, time_steps):
    trajectories = []
    for i in range(num_satellites):
        mean_anomaly = (2 * np.pi * i) / num_satellites  # Adjust mean anomaly for equidistant spacing
        trajectory = []
        for t in time_steps:
            true_anomaly = calculate_true_anomaly(semi_major_axis, eccentricity, mean_anomaly, t)
            x, y, z = satellite_position(semi_major_axis, eccentricity, np.deg2rad(inclination), np.deg2rad(argument_of_perigee), true_anomaly)
            
            # Calculate the gravitational force including the Moon's effect
            gravitational_force_total = gravitational_force_including_moon(masses[i], np.array([x, y, z]))
            
            # Calculate true anomaly with the Moon's gravitational effect
            true_anomaly += (gravitational_force_total / (masses[i] * semi_major_axis)) * t
            
            # Recalculate satellite position
            x, y, z = satellite_position(semi_major_axis, eccentricity, np.deg2rad(inclination), np.deg2rad(argument_of_perigee), true_anomaly)
            
            trajectory.append((x, y, z))
        trajectories.append(trajectory)
    #return np.array(trajectories)  # Convert to NumPy array
    return trajectories  # Convert to NumPy array

# Function to calculate true anomaly from time and orbital parameters
def calculate_true_anomaly(semi_major_axis, eccentricity, mean_anomaly, time):
    mean_motion = np.sqrt((G * (M_earth + M_moon)) / (semi_major_axis ** 3))
    mean_anomaly += mean_motion * time
    eccentric_anomaly = calculate_eccentric_anomaly(eccentricity, mean_anomaly)
    true_anomaly = 2 * np.arctan2(np.sqrt(1 + eccentricity) * np.sin(eccentric_anomaly / 2), np.sqrt(1 - eccentricity) * np.cos(eccentric_anomaly / 2))
    return true_anomaly

# Function to calculate eccentric anomaly using Newton-Raphson method
def calculate_eccentric_anomaly(eccentricity, mean_anomaly):
    tolerance = 1e-8
    max_iterations = 1000
    eccentric_anomaly = mean_anomaly
    for _ in range(max_iterations):
        next_eccentric_anomaly = eccentric_anomaly - (eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly) - mean_anomaly) / (1 - eccentricity * np.cos(eccentric_anomaly))
        if abs(next_eccentric_anomaly - eccentric_anomaly) < tolerance:
            return next_eccentric_anomaly
        eccentric_anomaly = next_eccentric_anomaly
    return eccentric_anomaly

# Satellite parameters
num_satellites = 8
semi_major_axis = 42164000.0  # Equidistant distance from Earth's center for geostationary orbit (approx. 42,164 km)
eccentricity = 0.1
inclination = 45  # Zero inclination for geostationary orbit (degrees)
argument_of_perigee = 0  # Zero argument of perigee for geostationary orbit (degrees)
masses = [100, 120, 90, 110, 130, 140, 150, 160]  # Masses of the satellites (kg)
time_steps = np.linspace(0, 3600, 100)  # Time steps in seconds (e.g., one hour)

# Calculate satellite trajectories including the Moon's gravitational effect
trajectories = calculate_satellite_trajectories_including_moon(num_satellites, semi_major_axis, eccentricity, inclination, argument_of_perigee, masses, time_steps)

'''
# Return the matrix of satellite trajectories
print(trajectories)

# Create a figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set axis labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Satellite Trajectories')

# Create a sphere to represent Earth
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = R_earth * np.outer(np.sin(u), np.cos(v))
y = R_earth * np.outer(np.sin(u), np.sin(v))
z = R_earth * np.outer(np.cos(u), np.ones(np.size(v)))

ax.plot_surface(x, y, z, color='b', alpha=0.1)

# Plot satellite trajectories
for i, trajectory in enumerate(trajectories):
    xs, ys, zs = zip(*trajectory)
    ax.plot(xs, ys, zs, label=f'Satellite {i+1}')

plt.legend()
plt.show()
'''