import numpy as np
import matplotlib.pyplot as plt

# Initialise Variables
T_initial = 60  # Initial geyser temperature in Celsius
T_cold = 25    # Cold water temperature in Celsius
T_ambient = 25  # Ambient temperature in Celsius (environment temp)
T_set = 60      # Desired maximum temperature
T_threshold = 55  # Threshold below which the heating element turns on
reheating_rate_per_minute = 2  # Heater raises temperature by 2°C per minute
reheating_rate_per_time_step = reheating_rate_per_minute * 0.6  # Adjusted for 0.01-hour time steps (0.6 minutes)

# Geyser Physical Properties
water_capacity = 200  # Capacity of geyser in liters
m = water_capacity  # Mass of the water in kg (1 liter = 1 kg for water)
cp = 4.186  # Specific heat capacity of water in J/kg°C

surface_area = 0.46

# Heat Transfer Coefficients
U = 1.2  # Overall heat transfer coefficient (in kJ/m^2°C)

k = 0.8  # Thermal conductivity of geyser walls - 1 (W/m·K)
dx = 0.009 # Thickness of the geyser walls in meters

# Time variables
time = np.arange(0, 24, 0.01)  # Time in hours (simulating a day with smaller time steps)
hot_water_usage = np.zeros_like(time)  # Empty array for hot water usage

# Simulate Hot Water Usage
for i, t in enumerate(time):
    if 7 <= t <= 8 or 19 <= t <= 20:
        hot_water_usage[i] = 40  # Use 20 liters during peak times

# Geyser Temperature over Time
temperature = np.full_like(time, T_initial)
heater_on = False  # Track whether the heater is on or off

# Function for convection (heating element on)

# Function for conduction (heat loss)
def conduction_step(T_fluid, dt):
    Q_cond = k * surface_area * (T_fluid - T_ambient) / dx  # Conduction heat loss
    dT_cond = (Q_cond / (m * cp)) * dt  # Temperature change due to conductio
    return T_fluid - dT_cond

def convection_step(T_fluid, dt):
    Q_conv = 4000 # Convection heat input
    dT_conv = (Q_conv / (m * cp)) * dt  # Temperature change due to convection
    return T_fluid + dT_conv

# Simulation Loop
for i in range(1, len(time)):
    # Calculate new temperature when hot water is used
    if hot_water_usage[i] > 0:
        used_water = hot_water_usage[i]
        inflow_water = used_water
        temperature[i] = ((temperature[i-1] * (water_capacity - inflow_water)) + (T_cold * inflow_water)) / water_capacity
    else:
        # Gradual heat loss through conduction when no water is used
        temperature[i] = conduction_step(temperature[i-1], 0.01)  # 0.01 hours per step

        # Heater turns on if temperature drops below the threshold
        if temperature[i] < T_threshold:
            heater_on = True

        # Heater is on: Reheat the water using convection
        if heater_on:
            temperature[i] = convection_step(temperature[i], 0.01)  # 0.01 hours per step
            if temperature[i] >= T_set:
                temperature[i] = T_set  # Ensure it caps at 60°C
                heater_on = False  # Turn off the heater
        
    print(temperature[i])
# Plot the Results
plt.figure(figsize=(12, 6))

# Hot Water Usage Plot
plt.subplot(2, 1, 1)
plt.plot(time, hot_water_usage, label='Hot Water Usage (litres)')
plt.ylabel('Hot Water Usage (litres)')
plt.xlabel('Time (hours)')
plt.legend()

# Geyser Temperature Plot
plt.subplot(2, 1, 2)
plt.plot(time, temperature, label='Geyser Temperature (°C)', color='r')
plt.ylabel('Temperature (°C)')
plt.xlabel('Time (hours)')
plt.legend()

plt.tight_layout()
plt.show()
