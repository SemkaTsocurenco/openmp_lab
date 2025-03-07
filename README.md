# openmp_lab

## Задание
Для выполнения работы требуется решить нестационарное уравнение теплопроводности для прямоугольной пластины, используя метод неявных разностных схем. Исходные условия для задачи включают следующие требования:


* **Начальная температура пластины:** равномерно распределена и составляет `T_init = 20 C`.
* **Граничные условия:**
    + Температура на нижней границе пластины фиксирована: `T_base = 50 C`.
    + Температура на наклонных границах поддерживается на уровне `T_side = 100 C`.
    
* **Геометрия пластины:**
   + Ширина пластины составляет `L_x = 8` единиц.
   + Высота пластины составляет `L_y = 4` единицы.
   + Треугольная область внутри прямоугольника определяется условием:  `y = -x + 8`  и `y = x` .


![изображение](https://github.com/user-attachments/assets/bcd9bafe-140a-4899-a491-ab3ce7c43dcf)

Основная цель заключается в реализации программы, моделирующей динамику температурного поля внутри указанной области. Для решения задачи потребуется:


* Построить численную схему, основанную на неявных разностных методах.
* Организовать эффективное вычисление температуры с использованием метода разделения переменных и распараллеливания вычислений на уровне обработки сетки (например, с помощью OpenMP).
* Реализовать функции для сохранения результатов температурного распределения для последующей визуализации с помощью Python-библиотеки matplotlib.


Дополнительно, нужно рассмотреть аспекты производительности:
  * Ускорение программы в зависимости от количества потоков.
  * Эффективность распараллеливания.
  * Сложность задачи и изоэффективность.

Графическое представление результатов включает визуализацию изменения температуры на треугольной области в разные временные моменты (пример представления показан на изображении).

![collage](https://github.com/user-attachments/assets/db231bc3-c262-4b5a-a184-f4c76a2e91d9)

Так же программа строит графики на основе данных полученных из модели (пример на изображении): 

![graf_2](https://github.com/user-attachments/assets/aed9dbdd-ba4c-4793-a241-ef1809c2d55c) ![graf_1](https://github.com/user-attachments/assets/7848dfb8-93d6-4f02-a928-fca2ccad63c6)

# Запуск программы
нужно запустить скрипт run.sh

> `sudo chmod +x ./run.sh && ./run.sh`



