#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

const int Nx = 100;     // Количество точек по оси X
const int Ny = 100;     // Количество точек по оси Y
const double Lx = 8.0; // Длина основания треугольника
const double Ly = 4.0; // Высота треугольника
const double alpha = 0.05; // Коэффициент теплопроводности
const double dt = 0.01;    // Шаг по времени
const double T_init = 20.0; // Начальная температура
const double T_boundary_side = 100.0; // Температура на наклонных сторонах
const double T_boundary_base = 50.0;  // Температура на основании
const int steps = 3000;      // Количество временных шагов
const int save_interval = 100; // Интервал сохранения данных для анимации

// Проверка, находится ли точка внутри равнобедренного треугольника
bool isInsideTriangle(double x, double y) {
    return (y >= 0) && (y <= -x + 8) && (y <= x);
}

// Функция для сохранения текущего распределения температуры в файл
void saveTemperature(const std::vector<std::vector<double>>& T, int step) {
    std::ofstream outFile("temperature_step_" + std::to_string(step) + ".dat");
    if (!outFile) {
        std::cerr << "Ошибка при открытии файла для записи." << std::endl;
        return;
    }
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double x = i * (Lx / (Nx - 1));
            double y = j * (Ly / (Ny - 1));
            if (isInsideTriangle(x, y)) {
                outFile << x << " " << y << " " << T[i][j] << "\n";
            }
        }
    }
    outFile.close();
}

int main() {
    std::vector<std::vector<double>> T(Nx, std::vector<double>(Ny, T_init));
    std::vector<std::vector<double>> T_new(Nx, std::vector<double>(Ny, T_init));
    double dx = Lx / (Nx - 1);
    double dy = Ly / (Ny - 1);

    // Основной цикл по времени
    for (int step = 0; step < steps; ++step) {
        // Расчет температуры внутри треугольника
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double x = i * dx;
                double y = j * dy;

                if (isInsideTriangle(x, y)) {
                    double dTdx2 = (T[i + 1][j] - 2 * T[i][j] + T[i - 1][j]) / (dx * dx);
                    double dTdy2 = (T[i][j + 1] - 2 * T[i][j] + T[i][j - 1]) / (dy * dy);
                    T_new[i][j] = T[i][j] + alpha * dt * (dTdx2 + dTdy2);
                }
            }
        }

        // Явная фиксация температуры для граничных условий
        for (int i = 0; i < Nx; ++i) {
            double x = i * dx;
            double y_left = x;
            double y_right = -x + 8;

            // Устанавливаем температуру на наклонных сторонах
            if (y_left >= 0 && y_left < Ly) {
                T_new[i][int(y_left / dy)] = T_boundary_side;
                T[i][int(y_left / dy)] = T_boundary_side; // Фиксируем для T и T_new
            }
            if (y_right >= 0 && y_right < Ly) {
                T_new[i][int(y_right / dy)] = T_boundary_side;
                T[i][int(y_right / dy)] = T_boundary_side; // Фиксируем для T и T_new
            }
            T_new[i][1] = T_boundary_base;
            T[i][1] = T_boundary_base; // Фиксируем для T и T_new
        }


        // Обновляем массив температур
        std::swap(T, T_new);

        // Сохранение данных каждые `save_interval` шагов
        if (step % save_interval == 0) {
            saveTemperature(T, step);
            std::cout << "Сохранены данные для шага " << step << std::endl;
        }
    }

    return 0;
}
