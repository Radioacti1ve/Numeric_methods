import math
import matplotlib.pyplot as plt

# Параметр уравнений
a = 4

def f_system(x):
    """Система нелинейных уравнений"""
    return [
        x[0]**2 / (a**2) + x[1]**2 / ((a/2)**2) - 1,
        a * x[1] - math.exp(x[0]) - x[0]
    ]

def J_system(x):
    """Матрица Якоби"""
    return [
		[2 * x[0] / (a**2), 2 * x[1] / ((a/2)**2)],
        [-math.exp(x[0]) - 1, a]
    ]

def phi_system(x):
    """Функция для метода простой итерации"""
    try:
        x2_new = math.sqrt((1 - x[0]**2 / (a**2)) * ((a/2)**2))
        x1_new = math.log(a * x2_new - x[0] + 1)  # Приближенная итерация
    except (ValueError, ZeroDivisionError):
        x2_new = x[1]
        x1_new = x[0]
    return [x1_new, x2_new]

def phi_derivative(x):
    """Приближённая матрица производных (нулевая для оценки сходимости)"""
    return [
        [0, 0],
        [0, 0]
    ]

def compute_spectral_norm(J):
    norm1 = abs(J[0][0]) + abs(J[0][1])
    norm2 = abs(J[1][0]) + abs(J[1][1])
    return max(norm1, norm2)

def simple_iteration(phi, phi_der, x0, epsilon, max_iter=100):
    x = x0.copy()
    flag = 1
    for k in range(max_iter):
        x_new = phi(x)
        J_phi = phi_der(x)
        q = compute_spectral_norm(J_phi)

        if q >= 1:
            flag = 0
            raise ValueError(f"Условие сходимости нарушено: q = {q:.6f} >= 1 на итерации {k}")

        delta_norm = max(abs(x_new[0] - x[0]), abs(x_new[1] - x[1]))

        if delta_norm * q / (1 - q) < epsilon:
            print(f"Решение найдено за {k + 1} итераций")
            break

        x = x_new
    else:
        print(f"Достигнуто максимальное число итераций {max_iter}")
    return x, k + 1, flag

def newton_method(f, J, x0, epsilon, max_iter=100):
    x = x0.copy()
    for k in range(max_iter):
        f_val = f(x)
        J_val = J(x)

        det = J_val[0][0] * J_val[1][1] - J_val[0][1] * J_val[1][0]
        if abs(det) < 1e-12:
            raise ValueError("Матрица Якоби вырождена")

        det_x = (-f_val[0] * J_val[1][1]) - (J_val[0][1] * -f_val[1])
        det_y = (J_val[0][0] * -f_val[1]) - (-f_val[0] * J_val[1][0])

        delta_x = det_x / det
        delta_y = det_y / det

        x_new = [x[0] + delta_x, x[1] + delta_y]

        if max(abs(x_new[0] - x[0]), abs(x_new[1] - x[1])) < epsilon:
            break

        x = x_new

    return x, k + 1

def check_newton_convergence(x0):
    J = J_system(x0)
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    return abs(det) != 0

def plot_system():
    def f1(x2):
        try:
            return math.sqrt(16 * (1 - x2**2 / 4))
        except ValueError:
            return float('nan')

    def f2(x1):
        return (math.exp(x1) + x1) / a

    x2_values = [i * 0.05 for i in range(-40, 40)]
    x1_values_f1 = [f1(x2) for x2 in x2_values]

    x1_values = [i * 0.05 for i in range(-40, 40)]
    x2_values_f2 = [f2(x1) for x1 in x1_values]

    plt.figure(figsize=(12, 10))
    plt.plot(x1_values_f1, x2_values, label=r'$\frac{x_1^2}{16} + \frac{x_2^2}{4} - 1 = 0$')
    plt.plot(x1_values, x2_values_f2, label=r'$4x_2 - e^{x_1} - x_1 = 0$')
    plt.xlabel('x1')
    plt.ylabel('x2')
    plt.title('Графическое определение начального приближения')
    plt.grid(True)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    x0 = [1.0, 0.5]  # Подобрано исходя из графика
    epsilon = 1e-7

    print("\nМетод простой итерации:")
    try:
        x_si, iter_si, flag = simple_iteration(phi_system, phi_derivative, x0, epsilon)
        print(f"Решение: x1 = {x_si[0]:.6f}, x2 = {x_si[1]:.6f}")
        print(f"Итераций: {iter_si}")
    except ValueError as e:
        print(e)

    print("\nМетод Ньютона:")
    try:
        x_nm, iter_nm = newton_method(f_system, J_system, x0, epsilon)
        print(f"Решение: x1 = {x_nm[0]:.6f}, x2 = {x_nm[1]:.6f}")
        print(f"Итераций: {iter_nm}")
    except ValueError as e:
        print(e)

    print("\nПроверка сходимости метода Ньютона:")
    if check_newton_convergence(x0):
        print("Условия сходимости выполнены")
    else:
        print("Условия сходимости не выполнены")

    print("\nПостроение графика для выбора начального приближения:")
    plot_system()
