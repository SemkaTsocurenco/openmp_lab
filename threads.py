import pandas as pd
import matplotlib.pyplot as plt

# Чтение данных из файла
data = pd.read_csv("results.txt", sep=" ", names=["threads", "size", "time", "complexity"])

# Группировка данных по количеству потоков и построение графика ускорения
for size in sorted(data["size"].unique()):
    subset = data[data["size"] == size]
    T1 = subset[subset["threads"] == 1]["time"].values[0]
    subset["speedup"] = T1 / subset["time"]
    plt.plot(subset["threads"], subset["speedup"], marker='o', label=f"Size: {size}")

plt.xlabel("Число потоков")
plt.ylabel("Ускорение")
plt.title("Ускорение работы программы в зависимости от числа потоков")
plt.legend()
plt.grid()
plt.show()

# Группировка данных по сложности задачи и построение графика ускорения
for threads in sorted(data["threads"].unique()):
    subset = data[data["threads"] == threads]
    T1 = subset[subset["threads"] == 1]["time"].values[0]
    subset["speedup"] = T1 / subset["time"]
    plt.plot(subset["complexity"], subset["speedup"], marker='o', label=f"Threads: {threads}")

plt.xlabel("Сложность задачи")
plt.ylabel("Ускорение")
plt.title("Ускорение работы программы в зависимости от сложности задачи")
plt.legend()
plt.grid()
plt.show()
