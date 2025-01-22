import pandas as pd
import matplotlib.pyplot as plt

# Read data
data = pd.read_csv("results.txt", sep=" ", names=["threads", "size", "time", "complexity"])

# Ensure no "SettingWithCopyWarning"
data = data.copy()

# Plot speedup vs threads for each size
for size in sorted(data["size"].unique()):
    subset = data[data["size"] == size]
    T1_data = subset[subset["threads"] == 1]
    if not T1_data.empty:
        T1 = T1_data["time"].values[0]
        subset["speedup"] = T1 / subset["time"]
        plt.plot(subset["threads"], subset["speedup"], marker='o', label=f"Size: {size}")
    else:
        print(f"No data for threads=1 with size={size}, skipping.")

plt.xlabel("Число потоков")
plt.ylabel("Ускорение")
plt.title("Ускорение работы программы в зависимости от числа потоков")
plt.legend()
plt.grid()
plt.savefig ("./graf_1.jpg")
plt.show()

# Plot speedup vs complexity for each number of threads
for threads in sorted(data["threads"].unique()):
    subset = data[data["threads"] == threads]
    T1_data = data[data["threads"] == 1]  # Use global T1 for consistency
    if not T1_data.empty:
        T1 = T1_data["time"].values[0]
        subset["speedup"] = T1 / subset["time"]
        plt.plot(subset["complexity"], subset["speedup"], marker='o', label=f"Threads: {threads}")
    else:
        print(f"No data for threads=1, skipping.")

plt.xlabel("Сложность задачи")
plt.ylabel("Ускорение")
plt.title("Ускорение работы программы в зависимости от сложности задачи")
plt.legend()
plt.grid()
plt.savefig ("./graf_2.jpg")
plt.show()

