def get_index(x, x0):
    for i in range(len(x)):
        if abs(x[i] - x0) < 1e-8:
            return i
    return -1

def first_derivative(x, y, index_x):
    if 0 < index_x < len(x) - 1:
        print("Точка находится внутри отрезка: ")
        dx_left = (y[index_x] - y[index_x - 1]) / (x[index_x] - x[index_x - 1])
        dx_right = (y[index_x + 1] - y[index_x]) / (x[index_x + 1] - x[index_x])
        print(f"  Левосторонняя производная: {dx_left:.6f}")
        print(f"  Правосторонняя производная: {dx_right:.6f}")
        diff_dy = (dx_right - dx_left) / (x[index_x + 1] - x[index_x - 1])
        derivative = dx_left + diff_dy * (2 * x[index_x] - x[index_x - 1] - x[index_x])
        print(f"  Первая производная: {derivative:.6f} \n")
    elif index_x == 0:
        print("Точка находится на левом крае:")
        dx_right = (y[index_x + 1] - y[index_x]) / (x[index_x + 1] - x[index_x])
        print(f"  Правая производная: {dx_right:.6f} \n")
    elif index_x == len(x) - 1:
        print("Точка находится на правом крае: ")
        dx_left = (y[index_x] - y[index_x - 1]) / (x[index_x] - x[index_x - 1])
        print(f"  Левая производная: {dx_left:.6f} \n")

def second_derivative_standard(x, y, index_x):
    if 0 < index_x < len(x) - 1:
        h = x[index_x + 1] - x[index_x]  # шаг (у нас он равен 0.2)
        dx2 = (y[index_x - 1] - 2 * y[index_x] + y[index_x + 1]) / (h ** 2)
        print(f"Вторая производная (центральная): {dx2:.6f}")
    else:
        print("Вторая производная не определена на краях.")

def second_derivative(x, y, index_x):
    if 0 < index_x < len(x) - 1:
        dx2 = 2 * ((y[index_x + 1] - y[index_x]) / (x[index_x + 1] - x[index_x]) -
                   (y[index_x] - y[index_x - 1]) / (x[index_x] - x[index_x - 1])) / \
              (x[index_x + 1] - x[index_x - 1])
        print(f"Вторая производная: {dx2:.6f}")
    else:
        print("Вторая производная не определена на краях.")

# --- данные из таблицы варианта 15 ---
x  = [0.0, 0.2, 0.4, 0.6, 0.8]
y  = [1.0, 1.4214, 1.8918, 2.4221, 3.0255]
x0 = 0.4  # X*

index_x = get_index(x, x0)
if index_x != -1:
    first_derivative(x, y, index_x)
    second_derivative_standard(x, y, index_x)
    second_derivative(x, y, index_x)
else:
    print("Точка не найдена в таблице")
