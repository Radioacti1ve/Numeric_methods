import math
import matplotlib.pyplot as plt
from typing import Callable

# равномерная сетка
def linspace(a: float, b: float, h: float) -> list[float]:
    n = int(round((b - a) / h)) + 1
    return [a + i * h for i in range(n)]

# приводим ОДУ 2-го порядка к системе:
# y' = z
# из x*z' + z = 0  ->  z' = -z/x
def f(x: float, y: float, z: float) -> float:
    return z

def g(x: float, y: float, z: float) -> float:
    # деление на 0 на нашем отрезке не встречается (x ∈ [1,2]),
    # оставим защиту на случай изменения интервала
    if x == 0:
        raise ZeroDivisionError("x=0 недопустим для уравнения z' = -z/x")
    return -z / x

# точное решение: y = 1 + ln(x)
def y_exact(x: float) -> float:
    return 1.0 + math.log(x)

class ODE2_Euler:
    def __init__(self, x0, xn, h, y0, z0, f, g):
        self.X = linspace(x0, xn, h)
        n = len(self.X)
        self.Y = [y0] * n
        self.Z = [z0] * n
        for i in range(1, n):
            self.Z[i] = self.Z[i-1] + h * g(self.X[i-1], self.Y[i-1], self.Z[i-1])
            self.Y[i] = self.Y[i-1] + h * f(self.X[i-1], self.Y[i-1], self.Z[i-1])

    def MAE(self, true: Callable):
        return max(abs(self.Y[i] - true(self.X[i])) for i in range(len(self.X)))

class ODE2_Runge_Kutta:
    def __init__(self, x0, xn, h, y0, z0, f, g):
        self.X = linspace(x0, xn, h)
        n = len(self.X)
        self.Y = [y0] * n
        self.Z = [z0] * n
        for i in range(n - 1):
            xi, yi, zi = self.X[i], self.Y[i], self.Z[i]
            K1 = h * f(xi, yi, zi)
            L1 = h * g(xi, yi, zi)
            K2 = h * f(xi + h/2, yi + K1/2, zi + L1/2)
            L2 = h * g(xi + h/2, yi + K1/2, zi + L1/2)
            K3 = h * f(xi + h/2, yi + K2/2, zi + L2/2)
            L3 = h * g(xi + h/2, yi + K2/2, zi + L2/2)
            K4 = h * f(xi + h, yi + K3, zi + L3)
            L4 = h * g(xi + h, yi + K3, zi + L3)
            self.Y[i+1] = yi + (K1 + 2*K2 + 2*K3 + K4) / 6
            self.Z[i+1] = zi + (L1 + 2*L2 + 2*L3 + L4) / 6

    def MAE(self, true: Callable):
        return max(abs(self.Y[i] - true(self.X[i])) for i in range(len(self.X)))

class ODE2_Adams:
    def __init__(self, x0, xn, h, y0, z0, f, g):
        self.X = linspace(x0, xn, h)
        n = len(self.X)
        if n < 4:
            raise ValueError("Для метода Адамса нужно минимум 4 точки")
        self.Y = [y0] * n
        self.Z = [z0] * n
        F = [0.0] * n
        G = [0.0] * n
        # стартовые 4 точки по Рунге–Кутте
        rk = ODE2_Runge_Kutta(x0, x0 + 3*h, h, y0, z0, f, g)
        self.Y[:4] = rk.Y[:4]
        self.Z[:4] = rk.Z[:4]
        for i in range(4):
            F[i] = f(self.X[i], self.Y[i], self.Z[i])
            G[i] = g(self.X[i], self.Y[i], self.Z[i])
        for i in range(3, n - 1):
            self.Y[i+1] = self.Y[i] + h/24 * (55*F[i] - 59*F[i-1] + 37*F[i-2] - 9*F[i-3])
            self.Z[i+1] = self.Z[i] + h/24 * (55*G[i] - 59*G[i-1] + 37*G[i-2] - 9*G[i-3])
            F[i+1] = f(self.X[i+1], self.Y[i+1], self.Z[i+1])
            G[i+1] = g(self.X[i+1], self.Y[i+1], self.Z[i+1])

    def MAE(self, true: Callable):
        return max(abs(self.Y[i] - true(self.X[i])) for i in range(len(self.X)))

class RungeRomberg:
    def __init__(self, X, Y_h, Y_h2, p):
        self.X = X
        self.Y = [y2 + (y2 - y1) / (2**p - 1) for y1, y2 in zip(Y_h, Y_h2[::2])]

    def MAE(self, true: Callable):
        return max(abs(self.Y[i] - true(self.X[i])) for i in range(len(self.X)))

def plot(solvers, labels, title):
    plt.figure(figsize=(10, 6))
    x_exact = linspace(1, 2, 0.01)
    plt.plot(x_exact, [y_exact(xi) for xi in x_exact], label="Точное решение", linewidth=2)
    for solver, label in zip(solvers, labels):
        plt.plot(solver.X, solver.Y, '--', label=label)
    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.show()

# Параметры вашей задачи
x0, xn, h = 1.0, 2.0, 0.1
y0, z0 = 1.0, 1.0  # y(1)=1, y'(1)=1

# Эйлер
euler1 = ODE2_Euler(x0, xn, h, y0, z0, f, g)
euler2 = ODE2_Euler(x0, xn, h/2, y0, z0, f, g)
rr_euler = RungeRomberg(euler1.X, euler1.Y, euler2.Y, 2)

# Рунге–Кутта 4-го порядка
runge1 = ODE2_Runge_Kutta(x0, xn, h, y0, z0, f, g)
runge2 = ODE2_Runge_Kutta(x0, xn, h/2, y0, z0, f, g)
rr_runge = RungeRomberg(runge1.X, runge1.Y, runge2.Y, 4)

# Адамс–Башфорт (4 шага)
adams1 = ODE2_Adams(x0, xn, h, y0, z0, f, g)
adams2 = ODE2_Adams(x0, xn, h/2, y0, z0, f, g)
rr_adams = RungeRomberg(adams1.X, adams1.Y, adams2.Y, 4)

# Сборка и подписи
solvers = [s for s in [euler1, runge1, adams1] if s is not None]
labels  = [l for l, s in zip(["Эйлер", "Рунге–Кутта", "Адамс"], solvers) if s is not None]

# Ошибки
print("Метод Эйлера:")
print(f"MAE: {euler1.MAE(y_exact):.10f}")
print(f"Runge–Romberg MAE: {rr_euler.MAE(y_exact):.10f}\n")

print("Метод Рунге–Кутта:")
print(f"MAE: {runge1.MAE(y_exact):.10f}")
print(f"Runge–Romberg MAE: {rr_runge.MAE(y_exact):.10f}\n")

print("Метод Адамса:")
print(f"MAE: {adams1.MAE(y_exact):.10f}")
print(f"Runge–Romberg MAE: {rr_adams.MAE(y_exact):.10f}\n")

# Графики
plot([euler1, rr_euler], ["Эйлер", "Эйлер (Рунге–Ромберг)"], "Метод Эйлера")
plot([runge1, rr_runge], ["Рунге–Кутта", "Рунге–Кутта (Рунге–Ромберг)"], "Метод Рунге–Кутты")
plot([adams1, rr_adams], ["Адамс", "Адамс (Рунге–Ромберг)"], "Метод Адамса")
plot(solvers, labels, "Сравнение методов")