import numpy as np
from scipy.integrate import odeint

# Define parameters
k = 0.2  # Thermal conductivity of PVC (W/m*K)
rho = 1.38  # Density (kg/m^3)
cp = 1339  # Specific heat capacity (J/kg*K)
r_inner = 0.03  # Inner radius of the oven (m)
r_outer = 0.04  # Outer radius of the oven (m)
h = 10  # Heat transfer coefficient at the outer surface (W/m^2*K)
T_inf = 25  # Ambient temperature (°C)
T_initial = 150  # Initial internal temperature (°C)

# Define the time span (2 hours = 7200 seconds)
t_span = np.linspace(0, 7200, 1000)

# Define the heat conduction equation
def heat_conduction(T, t):
    # Calculate the temperature gradient in the radial direction
    dT_dr = np.gradient(T, r)
    
    # Calculate the heat transfer at the outer surface
    q = -k * 2 * np.pi * r_outer * h * (T - T_inf)
    
    # Calculate the change in temperature with time
    dT_dt = (1 / (rho * cp)) * (1 / r) * np.gradient(r * dT_dr, r) + q / (rho * cp)
    
    return dT_dt

# Create an array of radial positions (discretization)
r = np.linspace(r_inner, r_outer, 100)

# Set initial temperature distribution (constant initial temperature)
T0 = np.full(len(r), T_initial)

# Solve the heat conduction equation over time
T_solution = odeint(heat_conduction, T0, t_span)

# Calculate the energy usage (integrate heat transfer over time)
# This assumes that the energy usage is directly proportional to the heat transfer
energy_used = np.trapz(np.trapz(T_solution - T_initial, r), t_span)  # in Joules

# Convert energy to kilowatt-hours (kWh)
energy_used_kwh = energy_used / (3600 * 1000)

print(f"Energy used to maintain the internal temperature: {energy_used_kwh:.4f} kWh")
