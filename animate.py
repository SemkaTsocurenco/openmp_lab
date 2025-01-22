import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from matplotlib.ticker import MaxNLocator

def read_data(filename):
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    return x, y, z

# Сортировка файлов по числовой части имени
files = sorted(glob.glob("./potok1_size90/temperature_step_*.dat"), key=lambda x: int(re.findall(r'\d+', x)[-1]))

# Определяем количество строк и столбцов в коллаже
rows = 5
cols = len(files) // rows + (len(files) % rows > 0)

# Создаём фигуру для коллажа
fig, axes = plt.subplots(rows, cols, figsize=(30, 30))
axes = axes.flatten()

for i, file in enumerate(files):
    x, y, z = read_data(file)
    ax = axes[i]
    contour = ax.tricontourf(x, y, z, levels=100, cmap='hot')
    ax.set_title(f"Шаг времени: {i * 100}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Добавляем цветовую шкалу
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Temperature (°C)')

# Отключаем пустые оси, если файлов меньше, чем ячеек в коллаже
for j in range(len(files), len(axes)):
    fig.delaxes(axes[j])

# Устанавливаем общий заголовок
fig.suptitle("Коллаж распределения температуры на пластине", fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig("./collage.jpg")
plt.show()

