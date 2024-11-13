import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from matplotlib.animation import FuncAnimation

def read_data(filename):
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    return x, y, z

files = sorted(glob.glob("./potok1_size90/temperature_step_*.dat"), key=lambda x: int(re.findall(r'\d+', x)[-1]))

fig, ax = plt.subplots()
contour = None
cbar = None

def update(frame):
    global contour, cbar
    if contour:
        for coll in contour.collections:
            coll.remove()
    x, y, z = read_data(files[frame])
    contour = ax.tricontourf(x, y, z, levels=100, cmap='hot')

    if cbar is None:
        cbar = plt.colorbar(contour, ax=ax, label='Temperature (°C)')
    ax.set_title(f"Шаг времени: {frame * 100}")
    return contour

ani = FuncAnimation(fig, update, frames=len(files), repeat=True)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Анимация распределения температуры на треугольной пластине')
plt.show()
