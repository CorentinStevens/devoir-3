import matplotlib.pyplot as plt
import numpy as np

# Read the data from the file
data = np.loadtxt('./time.txt')

# Extract columns
time = data[:, 0]
ux = data[:, 1]
uy = data[:, 2]
vx = data[:, 3]
vy = data[:, 4]

# Plot displacement (ux, uy) vs time
plt.figure(figsize=(10, 6))
plt.plot(time, ux, label='ux (Displacement X)', color='blue')
plt.plot(time, uy, label='uy (Displacement Y)', color='green')
plt.xlabel('Time (s)')
plt.ylabel('Displacement')
plt.title('Displacement vs Time')
plt.legend()
plt.grid()
plt.show()

# Plot velocity (vx, vy) vs time
plt.figure(figsize=(10, 6))
plt.plot(time, vx, label='vx (Velocity X)', color='red')
plt.plot(time, vy, label='vy (Velocity Y)', color='orange')
plt.xlabel('Time (s)')
plt.ylabel('Velocity')
plt.title('Velocity vs Time')
plt.legend()
plt.grid()
plt.show()