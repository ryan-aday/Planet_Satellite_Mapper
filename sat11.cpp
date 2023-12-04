#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

// Constants
const double G = 6.67430e-11;  // Gravitational constant (m^3/kg/s^2)
const double M_earth = 5.972e24;  // Earth's mass (kg)
const double M_moon = 7.342e22;   // Moon's mass (kg)
const double R_earth = 6.371e6;  // Earth's radius (m)
const double R_moon = 1.737e6;   // Moon's radius (m)

// Structure for 3D vectors
struct Vector3D {
    double x, y, z;
};

// Function to calculate the gravitational force on a satellite including the Moon's effect
Vector3D gravitationalForceIncludingMoon(double massSatellite, Vector3D r) {
    // Calculate the gravitational forces from both the Earth and the Moon
    double forceEarth = (G * M_earth * massSatellite) / (r.x * r.x + r.y * r.y + r.z * r.z);
    Vector3D moonPosition = {0, 0, 384400e3};  // Moon's position relative to Earth
    Vector3D moonToSatellite = {r.x - moonPosition.x, r.y - moonPosition.y, r.z - moonPosition.z};
    double distanceToMoon = sqrt(moonToSatellite.x * moonToSatellite.x + moonToSatellite.y * moonToSatellite.y + moonToSatellite.z * moonToSatellite.z);
    double forceMoon = (G * M_moon * massSatellite) / (distanceToMoon * distanceToMoon);

    // Calculate the direction of the gravitational force from the Moon
    Vector3D directionToMoon = {moonToSatellite.x / distanceToMoon, moonToSatellite.y / distanceToMoon, moonToSatellite.z / distanceToMoon};

    // Calculate the net gravitational force including both the Earth and the Moon
    Vector3D gravitationalForceTotal;
    gravitationalForceTotal.x = forceEarth * r.x / sqrt(r.x * r.x + r.y * r.y + r.z * r.z) - forceMoon * directionToMoon.x;
    gravitationalForceTotal.y = forceEarth * r.y / sqrt(r.x * r.x + r.y * r.y + r.z * r.z) - forceMoon * directionToMoon.y;
    gravitationalForceTotal.z = forceEarth * r.z / sqrt(r.x * r.x + r.y * r.y + r.z * r.z) - forceMoon * directionToMoon.z;

    return gravitationalForceTotal;
}

// Function to convert Cartesian coordinates to longitude, latitude, and altitude
void cartesianToGeodetic(double x, double y, double z, double& longitude, double& latitude, double& altitude) {
    // WGS84 parameters
    const double a = 6378137.0;  // Semi-major axis (m)
    const double b = 6356752.3142;  // Semi-minor axis (m)
    
    // Calculate other parameters
    double p = sqrt(x * x + y * y);
    double e2 = 1 - (b * b) / (a * a);
    
    double ep = sqrt((a * a - b * b) / (b * b));
    double th = atan2(a * z, b * p);
    
    // Calculate longitude
    longitude = atan2(y, x);

    // Calculate latitude
    double latA = z + ep * ep * b * pow(sin(th), 3);
    double latB = p - e2 * a * pow(cos(th), 3);
    latitude = atan2(latA, latB);
    
    // Calculate altitude
    double N = a / sqrt(1 - e2 * pow(sin(latitude), 2));
    altitude = p / cos(latitude) - N;

    // Convert latitude and longitude from radians to degrees
    latitude = latitude * 180.0 / M_PI;
    longitude = longitude * 180.0 / M_PI;
}

// Function to calculate the satellite's position in orbit
Vector3D satellitePosition(double semiMajorAxis, double eccentricity, double inclination, double argumentOfPerigee, double trueAnomaly) {
    double r = semiMajorAxis * (1 - eccentricity * eccentricity) / (1 + eccentricity * cos(trueAnomaly));
    Vector3D position;
    position.x = r * (cos(argumentOfPerigee) * cos(trueAnomaly) - sin(argumentOfPerigee) * sin(trueAnomaly) * cos(inclination));
    position.y = r * (sin(argumentOfPerigee) * cos(trueAnomaly) + cos(argumentOfPerigee) * sin(trueAnomaly) * cos(inclination));
    position.z = r * (sin(trueAnomaly) * sin(inclination));
    return position;
}

// Function to calculate satellite trajectories over time including the Moon's gravitational effect
std::vector<std::vector<Vector3D>> calculateSatelliteTrajectoriesIncludingMoon(int numSatellites, double semiMajorAxis, double eccentricity, double inclination, double argumentOfPerigee, const std::vector<double>& masses, const std::vector<double>& timeSteps) {
    std::vector<std::vector<Vector3D>> trajectories;
    for (int i = 0; i < numSatellites; ++i) {
        double meanAnomaly = (2 * M_PI * i) / numSatellites;  // Adjust mean anomaly for equidistant spacing
        std::vector<Vector3D> trajectory;
        for (double t : timeSteps) {
            double trueAnomaly = meanAnomaly + sqrt((G * (M_earth + M_moon)) / (semiMajorAxis * semiMajorAxis * semiMajorAxis)) * t;
            double r = (semiMajorAxis * (1 - eccentricity * eccentricity)) / (1 + eccentricity * cos(trueAnomaly));
            Vector3D position = satellitePosition(semiMajorAxis, eccentricity, inclination, argumentOfPerigee, trueAnomaly);
            
            // Calculate the gravitational force including the Moon's effect
            Vector3D gravitationalForceTotal = gravitationalForceIncludingMoon(masses[i], position);
            
            // Calculate true anomaly with the Moon's gravitational effect
            trueAnomaly += (sqrt((G * (M_earth + M_moon)) / (semiMajorAxis * semiMajorAxis * semiMajorAxis)) * t);
            
            // Recalculate satellite position
            position = satellitePosition(semiMajorAxis, eccentricity, inclination, argumentOfPerigee, trueAnomaly);
            
            trajectory.push_back(position);
        }
        trajectories.push_back(trajectory);
    }
    return trajectories;
}

int main() {
    // Satellite parameters
    int numSatellites = 8;
    double semiMajorAxis = 42164000.0;  // Equidistant distance from Earth's center for geostationary orbit (approx. 42,164 km)
    double eccentricity = 0.1;
    double inclination = 45.0;  // Zero inclination for geostationary orbit (degrees)
    double argumentOfPerigee = 0.0;  // Zero argument of perigee for geostationary orbit (degrees)
    std::vector<double> masses = {100, 120, 90, 110, 130, 140, 150, 160};  // Masses of the satellites (kg)

    std::vector<double> timeSteps;
    for (double t = 0.0; t <= 3600.0; t += 3600.0 / 100.0) {  // Time steps in seconds (e.g., one hour)
        timeSteps.push_back(t);
    }

    // Calculate satellite trajectories including the Moon's gravitational effect
    std::vector<std::vector<Vector3D>> trajectories = calculateSatelliteTrajectoriesIncludingMoon(numSatellites, semiMajorAxis, eccentricity, inclination, argumentOfPerigee, masses, timeSteps);

    // Output satellite trajectories to a text file
    std::ofstream outputFile("satellite_trajectories.txt");
    outputFile << std::fixed << std::setprecision(6);
    for (int i = 0; i < numSatellites; ++i) {
        for (const Vector3D& position : trajectories[i]) {
            outputFile << position.x << " " << position.y << " " << position.z << "\n";
        }
        outputFile << std::endl;
    }
    outputFile.close();

    // Output satellite trajectories in WGS84 coordinates to a text file
    std::ofstream wgs84OutputFile("satellite_trajectories_wgs84.txt");
    wgs84OutputFile << std::fixed << std::setprecision(6);
    for (int i = 0; i < numSatellites; ++i) {
        for (const Vector3D& position : trajectories[i]) {
            double longitude, latitude, altitude;
            cartesianToGeodetic(position.x, position.y, position.z, longitude, latitude, altitude);
            wgs84OutputFile << longitude << " " << latitude << " " << altitude << "\n";
        }
        wgs84OutputFile << std::endl;
    }
    wgs84OutputFile.close();

    return 0;
}
