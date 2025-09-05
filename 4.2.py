import math
import matplotlib.pyplot as plt

# ---- точное решение ----------------------------------------------------------
def exact_solution(x: float) -> float:
    # y(x) = 1 + x + ln x, определено только для x > 0
    if x <= 0:
        return float('nan')
    return 1.0 + x + math.log(x)

# ---- система первого порядка для RK4 ----------------------------------------
# y' = z
def f(x, y, z):
    return z

# z' = (x*z - y)/(x^2*ln x)
def g(x, y, z):
    if x <= 0:
        raise ValueError("x must be > 0 for ln(x).")
    ln = math.log(x)
    if abs(ln) < 1e-12:
        # Сингулярность в x=1; избегаем такого x выбором отрезка
        raise ZeroDivisionError("ln(x) ~ 0 → сингулярная точка (x≈1).")
    return (x * z - y) / (x * x * ln)

# ---- классический RK4 для системы y' = f, z' = g ----------------------------
def runge_kutta_4(f, g, a, b, y0, z0, h):
    n = int(round((b - a) / h))
    x = [a + i * h for i in range(n + 1)]
    y = [0.0] * (n + 1)
    z = [0.0] * (n + 1)
    y[0] = y0
    z[0] = z0

    for i in range(n):
        k1y = h * f(x[i], y[i], z[i])
        k1z = h * g(x[i], y[i], z[i])

        k2y = h * f(x[i] + h / 2, y[i] + k1y / 2, z[i] + k1z / 2)
        k2z = h * g(x[i] + h / 2, y[i] + k1y / 2, z[i] + k1z / 2)

        k3y = h * f(x[i] + h / 2, y[i] + k2y / 2, z[i] + k2z / 2)
        k3z = h * g(x[i] + h / 2, y[i] + k2y / 2, z[i] + k2z / 2)

        k4y = h * f(x[i] + h, y[i] + k3y, z[i] + k3z)
        k4z = h * g(x[i] + h, y[i] + k3y, z[i] + k3z)

        y[i + 1] = y[i] + (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0
        z[i + 1] = z[i] + (k1z + 2 * k2z + 2 * k3z + k4z) / 6.0

    return x, y

# ---- метод стрельбы для условий Дирихле y(a)=alpha, y(b)=beta ----------------
def shooting_method(f, g, a, b, alpha, beta, h, tol=1e-8, max_iter=100):
    # подбираем начальный z0 (т.е. y'(a)) секущей
    eta0 = 0.0
    eta1 = 1.0

    # первая траектория
    _, y0 = runge_kutta_4(f, g, a, b, alpha, eta0, h)
    phi0 = y0[-1] - beta

    for _ in range(max_iter):
        _, y1 = runge_kutta_4(f, g, a, b, alpha, eta1, h)
        phi1 = y1[-1] - beta
        if abs(phi1) < tol:
            return runge_kutta_4(f, g, a, b, alpha, eta1, h)

        # шаг секущей
        denom = (phi1 - phi0)
        if abs(denom) < 1e-14:
            denom = 1e-14
        eta_next = eta1 - phi1 * (eta1 - eta0) / denom
        eta0, eta1 = eta1, eta_next
        phi0 = phi1

    # если не сошлось — вернём последнюю траекторию
    return runge_kutta_4(f, g, a, b, alpha, eta1, h)

# ---- КР-схема (вторая точность) для a(x)y''+b(x)y'+c(x)y=d(x) ----------------
def finite_difference_method(aL, bR, alpha, beta, h):
    n = int(round((bR - aL) / h))
    x = [aL + i * h for i in range(n + 1)]

    def a_coef(xi):  # a(x) = x^2 ln x
        ln = math.log(xi)
        return xi * xi * ln

    def b_coef(xi):  # b(x) = -x
        return -xi

    def c_coef(xi):  # c(x) = 1
        return 1.0

    def d_rhs(xi):  # d(x) = 0
        return 0.0

    A = [0.0] * (n + 1)
    B = [0.0] * (n + 1)
    C = [0.0] * (n + 1)
    D = [0.0] * (n + 1)

    # ГУ Дирихле
    B[0] = 1.0
    D[0] = alpha
    B[-1] = 1.0
    D[-1] = beta

    # Внутренние узлы
    for i in range(1, n):
        xi = x[i]
        ln = math.log(xi)
        if abs(ln) < 1e-12:
            raise ZeroDivisionError("Попал в x≈1: ln(x) ~ 0. Выбери отрезок без x=1.")
        ai = a_coef(xi)
        bi = b_coef(xi)
        ci = c_coef(xi)
        di = d_rhs(xi)

        A[i] = ai / h**2 - bi / (2 * h)
        B[i] = -2 * ai / h**2 + ci
        C[i] = ai / h**2 + bi / (2 * h)
        D[i] = di

    # прогонка
    y = [0.0] * (n + 1)
    alpha_coef = [0.0] * (n + 1)
    beta_coef = [0.0] * (n + 1)

    alpha_coef[1] = -C[0] / B[0]
    beta_coef[1] = D[0] / B[0]

    for i in range(1, n):
        denom = B[i] + A[i] * alpha_coef[i]
        if abs(denom) < 1e-14:
            denom = 1e-14
        alpha_coef[i + 1] = -C[i] / denom
        beta_coef[i + 1] = (D[i] - A[i] * beta_coef[i]) / denom

    y[-1] = (D[-1] - A[-1] * beta_coef[-1]) / (B[-1] + A[-1] * alpha_coef[-1])

    for i in range(n - 1, -1, -1):
        y[i] = alpha_coef[i + 1] * y[i + 1] + beta_coef[i + 1]

    return x, y

# ---- Оценка Рунге—Ромберга ---------------------------------------------------
def compute_runge_romberg(x1, y1, x2, y2, p):
    # ищем общие узлы (при кратных шагах это совпадение каждого второго узла)
    idx2 = {round(val, 12): i for i, val in enumerate(x2)}
    Rs = []
    for i1, xi in enumerate(x1):
        key = round(xi, 12)
        if key in idx2:
            i2 = idx2[key]
            Rs.append((y2[i2] - y1[i1]) / (2 ** p - 1))
    if not Rs:
        return 0.0, []
    maxR = max(abs(r) for r in Rs)
    return maxR, Rs

# ================== ПАРАМЕТРЫ ЗАДАЧИ ==========================================
# Выбираем отрезок, не пересекающий x=1 (сингулярность a(x)=x^2 ln x).
a = 0.2
b = 0.9

# ГУ Дирихле из точного решения
alpha = exact_solution(a)
beta = exact_solution(b)

# Шаги сетки
h1 = 0.01
h2 = h1 / 2

# ================== ВЫЧИСЛЕНИЯ ===============================================

# Метод стрельбы (RK4)
x_rk_h1, y_rk_h1 = shooting_method(f, g, a, b, alpha, beta, h1)
x_rk_h2, y_rk_h2 = shooting_method(f, g, a, b, alpha, beta, h2)

# Конечно-разностный метод (2-й порядок)
x_fd_h1, y_fd_h1 = finite_difference_method(a, b, alpha, beta, h1)
x_fd_h2, y_fd_h2 = finite_difference_method(a, b, alpha, beta, h2)

# Точное решение на сетке h1
y_exact_h1 = [exact_solution(xi) for xi in x_rk_h1]
err_rk = [abs(yy - ye) for yy, ye in zip(y_rk_h1, y_exact_h1)]

# Точное на сетке КР
y_exact_fd_h1 = [exact_solution(xi) for xi in x_fd_h1]
err_fd = [abs(yy - ye) for yy, ye in zip(y_fd_h1, y_exact_fd_h1)]

# Оценка Рунге—Ромберга
max_R_rk, R_rk = compute_runge_romberg(x_rk_h1, y_rk_h1, x_rk_h2, y_rk_h2, p=4)
max_R_fd, R_fd = compute_runge_romberg(x_fd_h1, y_fd_h1, x_fd_h2, y_fd_h2, p=2)

# ================== ГРАФИКИ ===================================================
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(x_rk_h1, y_rk_h1, 'b-', label='Стрельба (RK4)')
plt.plot(x_fd_h1, y_fd_h1, 'g-', label='Конечно-разностный')
plt.plot(x_rk_h1, y_exact_h1, 'r--', label='Точное')
plt.title('Сравнение решений')
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 2)
plt.semilogy(x_rk_h1, err_rk, 'b-', label='Стрельба (ошибка)')
plt.semilogy(x_fd_h1, err_fd, 'g-', label='Кон.-разности (ошибка)')
plt.title('Погрешности (лог. шкала)')
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 3)
if R_rk:
    plt.plot(x_rk_h1[:len(R_rk)], R_rk, 'b-', label='РР для RK4 (p=4)')
if R_fd:
    plt.plot(x_fd_h1[:len(R_fd)], R_fd, 'g-', label='РР для КР (p=2)')
plt.title('Оценка Рунге—Ромберга')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

print("\nИтоговые результаты:")
print(f"Макс. погрешность стрельбы:        {max(err_rk):.6e}")
print(f"Макс. погрешность кон. разностей:  {max(err_fd):.6e}")
print(f"Оценка РР (стрельба):              {max_R_rk:.6e}")
print(f"Оценка РР (кон. разности):         {max_R_fd:.6e}")
