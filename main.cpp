#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <omp.h>
#include <boost/filesystem.hpp>

const double Lx = 8.0;
const double Ly = 4.0;
const double alpha = 0.05;
const double dt = 0.01;
const double T_init = 20.0;
const double T_boundary_side = 100.0;
const double T_boundary_base = 50.0;
const int steps = 3000;
const int save_interval = 100;

bool isInsideTriangle(double x, double y) {
    return (y >= 0) && (y <= -x + 8) && (y <= x);
}
void saveTemperature(const std::vector<std::vector<double>>& T, int step, int Nx, int Ny, int num_threads) {
    boost::filesystem::create_directory("potok"+std::to_string(num_threads) + "_size"+std::to_string(Nx));
    std::ofstream outFile("./potok"+std::to_string(num_threads) + "_size"+std::to_string(Nx)+"/"
    +"temperature_step_" + std::to_string(step) + ".dat");
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

void saveTemperature(const std::vector<std::vector<double>>& T, int step,  int Nx, int Ny) {
    boost::filesystem::create_directory("NoParallel");
    std::ofstream outFile("./NoParallel/temperature_step_" + std::to_string(step) + ".dat");
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

double run_simulation(int Nx, int Ny, int num_threads) {
    std::vector<std::vector<double>> T(Nx, std::vector<double>(Ny, T_init));
    std::vector<std::vector<double>> T_new(Nx, std::vector<double>(Ny, T_init));
    double dx = Lx / (Nx - 1);
    double dy = Ly / (Ny - 1);

    omp_set_num_threads(num_threads);
    double start = omp_get_wtime();

    #pragma omp parallel
    {
        for (int step = 0; step < steps; ++step) {
            #pragma omp for
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

            #pragma omp for
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



            #pragma omp barrier
            #pragma omp single
            {
                std::swap(T, T_new);
            }
            #pragma omp master
            if (step % save_interval == 0) {
                saveTemperature(T_new, step, Nx, Ny, num_threads);
            }
        }
    }

    double end = omp_get_wtime();
    return end - start;
}


double run_no_parallel(int Nx, int Ny) {
    std::vector<std::vector<double>> T(Nx, std::vector<double>(Ny, T_init));
    std::vector<std::vector<double>> T_new(Nx, std::vector<double>(Ny, T_init));
    double dx = Lx / (Nx - 1);
    double dy = Ly / (Ny - 1);
    double start = omp_get_wtime();
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
            saveTemperature(T, step, Nx, Ny);
        }
    }
    double end = omp_get_wtime();
    return end - start;
}


int main() {
    std::ofstream results("results.txt");
    if (!results) {
        std::cerr << "Ошибка при открытии файла results.txt" << std::endl;
        return 1;
    }
    double Tmin = 100.0;
    int min_num_trheads;
    int min_size;
    for (int num_threads = 1; num_threads <= 16; num_threads *= 2) {
        for (int size = 30; size <= 100; size += 30) {
            double time = run_simulation(size, size, num_threads);
            double complexity = size * size * steps;
            results << num_threads << " " << size << " " << time << " " << complexity << "\n";
            if (time < Tmin) {
                Tmin = time;
                min_num_trheads = num_threads;
                min_size = size;
            }
            std::cout << "Потоки: " << num_threads << ", Размер: " << size << ", Время: " 
                    << time << " сек, Сложность: " << complexity << " - DONE\n";
        }
    }

    std::cout<<"NO PARALLEL - DONE\n";

    double T_no_parallel = run_no_parallel(min_size, min_size);

    double speedup = T_no_parallel / Tmin;
    double efficiency = speedup / min_num_trheads;
    double complexity = min_size * min_size * steps;
    double isoefficiency = complexity / (min_num_trheads * min_num_trheads);

    std::cout << "\nВремя выполнения без разпараллеливания: " << T_no_parallel<< "\n";
    std::cout<< "Минимальное время выолнения с потоками " << Tmin<< "\n";
    std::cout << "Ускорение: " << speedup << "\n";
    std::cout << "Эффективность: " << efficiency << "\n";
    std::cout << "Сложность задачи: " << complexity << " операций\n";
    std::cout << "Значение функции изоэффективности: " << isoefficiency << "\n";


    results.close();
    return 0;
}
