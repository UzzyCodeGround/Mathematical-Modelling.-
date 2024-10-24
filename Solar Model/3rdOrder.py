import numpy as np
import matplotlib.pyplot as plt

# Initialise Variables
T_initial = 60  # Initial geyser temperature in Celsius
T_cold = 20     # Cold water temperature in Celsius
T_ambient = 25  # Ambient temperature in Celsius (environment temp)
T_set = 60      # Desired maximum temperature
T_threshold = 55  # Threshold below which the heating element turns on
reheating_rate_per_minute = 2  # Increase the reheating rate to 2°C per minute
reheating_rate_per_time_step = reheating_rate_per_minute * 0.6  # Adjusted for 0.01-hour time steps (0.6 minutes)

# Geyser Physical Properties
water_capacity = 100  # Capacity of geyser in liters
height = 1.5  # Height of the geyser in meters
radius = 0.3  # Radius of the geyser in meters
surface_area = 2 * np.pi * radius * height  # Surface area in square meters

# Heat Transfer Coefficient
U = 1.5  # Overall heat transfer coefficient (in kJ/m^2°C)

# Time variables
time = np.arange(0, 24, 0.01)  # Time in hours (simulating a day with smaller time steps)
hot_water_usage = np.zeros_like(time)  # Empty array for hot water usage

# Simulate Hot Water Usage
for i, t in enumerate(time):
    if 7 <= t <= 8 or 19 <= t <= 20:  
        hot_water_usage[i] = 20  # Use 20 liters during peak times

# Geyser Temperature over Time
temperature = np.full_like(time, T_initial)

# Simulation Loop
for i in range(1, len(time)):
    # Calculate new temperature when hot water is used
    if hot_water_usage[i] > 0:
        used_water = hot_water_usage[i]
        inflow_water = used_water
        temperature[i] = ((temperature[i-1] * (water_capacity - inflow_water)) + (T_cold * inflow_water)) / water_capacity
    else:
        # Gradual heat loss when no water is used
        heat_loss = U * surface_area * (temperature[i-1] - T_ambient) * 0.01  # Adjusted for 0.01-hour time steps
        temperature[i] = temperature[i-1] - heat_loss

        # If temperature drops below threshold, turn on heating element until it reaches T_set
        if temperature[i] < T_set:
            temperature[i] += reheating_rate_per_time_step  # Gradual reheating step-by-step

    # Cap temperature at the set point to prevent overheating
    if temperature[i] > T_set:
        temperature[i] = T_set

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