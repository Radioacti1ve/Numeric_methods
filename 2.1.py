import math
import matplotlib.pyplot as plt

# Функция f(x) = sin(x) - x^2 + 1
def f(x):
    return math.sin(x) - x**2 + 1

# Первая производная
def df(x):
    return math.cos(x) - 2 * x

# Вторая производная
def d2f(x):
    return -math.sin(x) - 2

# Итерационная функция phi(x) = sqrt(sin(x) + 1)
def phi(x):
    value = math.sin(x) + 1
    if value < 0:
        raise ValueError("Отрицательный аргумент под корнем")
    return math.sqrt(value)

# Производная phi(x)
def phi_derivative(x):
    value = math.sin(x) + 1
    if value <= 0:
        raise ValueError("Отрицательный аргумент под корнем")
    return 0.5 * math.cos(x) / math.sqrt(value)

# Метод простой итерации
def simple_iteration(phi, x0, epsilon, a, b, max_iter=1000):
    x_prev = x0
    iterations = 0
    for _ in range(max_iter):
        x_next = phi(x_prev)
        q = check_simple_iteration_conditions(phi, x_next, a, b)
        if abs(x_next - x_prev) * q / (1 - q) < epsilon:
            break
        x_prev = x_next
        iterations += 1
    return x_next, iterations

# Проверка условий на концах интервала для метода простой итерации
def check_simple_iteration(a, b):
    try:
        da = abs(phi_derivative(a))
        db = abs(phi_derivative(b))
        if da < 1 and db < 1:
            return True
        else:
            print(f"Производные на концах интервала не удовлетворяют условию сходимости: |phi'(a)|={da}, |phi'(b)|={db}")
            return False
    except ValueError:
        print("Ошибка при вычислении производной phi на концах интервала")
        return False

# Проверка условий сходимости внутри интервала
def check_simple_iteration_conditions(phi, x0, a, b, q_max=0.99):
    num_points = 100
    h = (b - a) / num_points
    max_derivative = 0
    for i in range(num_points + 1):
        x = a + i * h
        try:
            derivative = phi_derivative(x)
            current_abs_derivative = abs(derivative)
            if current_abs_derivative > max_derivative:
                max_derivative = current_abs_derivative
        except ValueError:
            continue
    if max_derivative >= q_max:
        print(f"Условие сходимости не выполнено (max |phi'| = {max_derivative:.4f})")
        return False
    return max_derivative

# Метод Ньютона
def newton_method(f, df, x0, epsilon=1e-6, max_iter=100):
    x_prev = x0
    iterations = 0
    for _ in range(max_iter):
        fx = f(x_prev)
        dfx = df(x_prev)
        if abs(dfx) < 1e-12:
            raise ValueError("Производная слишком близка к нулю!")
        x_next = x_prev - fx / dfx
        if abs(x_next - x_prev) < epsilon:
            break
        x_prev = x_next
        iterations += 1
    return x_next, iterations

# Проверка условий для метода Ньютона
def check_newton_conditions(f, d2f, a, b):
    if f(a) * f(b) >= 0:
        print("Начальное приближение вне отрезка [a, b]")
        return False, None
    if abs(d2f(a) * f(a)) < df(a) ** 2:
        return True, a
    elif abs(d2f(b) * f(b)) < df(b) ** 2:
        return True, b
    else:
        print("Не выполняется условие сходимости метода Ньютона")
        return False, None

# ===== Проверка на интервале [0.5, 1.5] =====
a, b = 0.5, 1.5
epsilon = 1e-8

print("\nПроверка интервала [0.5, 1.5]:")
if check_simple_iteration(a, b):
    x0_si = a if abs(phi_derivative(a)) < abs(phi_derivative(b)) else b
    root_si, iter_si = simple_iteration(phi, x0_si, epsilon, a, b)
    print(f"Метод простой итерации: корень = {root_si:.8f}, итераций = {iter_si}")

newton_ok, x0_nm = check_newton_conditions(f, d2f, a, b)
if newton_ok:
    root_nm, iter_nm = newton_method(f, df, x0_nm, epsilon)
    print(f"Метод Ньютона: корень = {root_nm:.8f}, итераций = {iter_nm}")

# ===== Проверка на интервале [-1.0, -0.5] =====
a, b = -1.0, -0.5

print("\nПроверка интервала [-1.0, -0.5]:")
if check_simple_iteration(a, b):
    x0_si = a if abs(phi_derivative(a)) < abs(phi_derivative(b)) else b
    root_si, iter_si = simple_iteration(phi, x0_si, epsilon, a, b)
    print(f"Метод простой итерации: корень = {root_si:.8f}, итераций = {iter_si}")

newton_ok, x0_nm = check_newton_conditions(f, d2f, a, b)
if newton_ok:
    root_nm, iter_nm = newton_method(f, df, x0_nm, epsilon)
    print(f"Метод Ньютона: корень = {root_nm:.8f}, итераций = {iter_nm}")

# ===== График функции =====
x_values = [i * 0.01 for i in range(-100, 200)]  # от -1 до 2
y_values = [f(x) for x in x_values]

plt.plot(x_values, y_values, label=r'$f(x) = \sin(x) - x^2 + 1$')
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('График функции f(x) = sin(x) - x^2 + 1')
plt.grid(True)
plt.legend()
plt.show()
