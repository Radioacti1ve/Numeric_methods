import math
import numpy as np
import matplotlib.pyplot as plt

# --- Начальные параметры задачи ---
# x ∈ [0, π], t ∈ [0, 5]
# h — шаг по x, τ — шаг по времени
# σ используется для выбора шага по устойчивости (σ ~ Courant number)
x_begin, x_end = 0.0, math.pi
t_begin, t_end = 0.0, 5.0
h = 0.01
sigma = 0.01
tau = math.sqrt(sigma) * h  # [Математика] τ = √σ * h (условие устойчивости)

# --- Аналитическое решение ---
# [Математика] Решение вида u(x,t) = 0.5 * e^{-x} * sin(x) * sin(2t)
def analytic_solution_fun(x, t):
    return 0.5 * math.exp(-x) * math.sin(x) * math.sin(2.0 * t)

# Начальные условия: u(x,0) и u_t(x,0)
def psi0(x):  # u(x,0)
    return 0.0

def psi1(x):  # u_t(x,0)
    return math.exp(-x) * math.sin(x)

# --- Граничные условия (Dirichlet) ---
# u(0,t) = 0, u(π,t) = 0
def phi0(t):
    return 0.0

def phi1(t):
    return 0.0

# --- Формирование сеток по x и t ---
def make_grids(x_begin, x_end, h, t_begin, t_end, tau):
    # возвращает массивы x и t для расчётов
    x = np.arange(x_begin, x_end + h*0.5, h)
    t = np.arange(t_begin, t_end + tau*0.5, tau)
    return x, t

x, t = make_grids(x_begin, x_end, h, t_begin, t_end, tau)

# --- Формирует аналитическое решение на сетке (для сравнения с численным) ---
def build_analytical_from_grids(x, t):
    U = np.zeros((len(t), len(x)))
    for ni in range(len(t)):
        for i in range(len(x)):
            U[ni, i] = analytic_solution_fun(x[i], t[ni])
    return x, t, U

# --- Метрики ошибки ---
def max_abs_error(A, B):
    assert A.shape == B.shape
    return np.max(np.abs(A - B))  # максимальная абсолютная ошибка

def mean_abs_error(A, B):
    assert A.shape == B.shape
    return np.mean(np.abs(A - B))  # средняя абсолютная ошибка

# --- Визуализация решений в фиксированный момент времени ---
def plot_at_time(solutions, x, t, time):
    idx = np.argmin(np.abs(t - time))
    plt.figure(figsize=(10,5))
    for name, S in solutions.items():
        plt.plot(x, S[idx, :], label=name)
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.legend()
    plt.grid()
    plt.show()

# --- Визуализация ошибки во времени ---
def plot_error_vs_time(solutions, analytical_name, t):
    plt.figure(figsize=(10,5))
    for name, S in solutions.items():
        if name == analytical_name:
            continue
        errs = [np.max(np.abs(S[k,:] - solutions[analytical_name][k,:])) for k in range(len(t))]
        plt.plot(t, errs, label=name)
    plt.xlabel('time')
    plt.ylabel('max abs error')
    plt.legend()
    plt.grid()
    plt.show()

# --- Явная конечно-разностная схема ---
# [Математика] u^{n+1}_i = 2u^n_i - u^{n-1}_i + τ²(2u_xx + 4u_x)
def explicit_scheme_from_grids(x, t, h, tau):
    nx = len(x); nt = len(t)
    U = np.zeros((nt, nx))

    # начальные условия
    for i in range(nx):
        U[0, i] = psi0(x[i])
        U[1, i] = psi0(x[i]) + tau * psi1(x[i])  # аппроксимация u(x,τ)

    # основной цикл по времени
    for n in range(1, nt-1):
        # граничные условия
        U[n+1, 0] = phi0(t[n+1])
        U[n+1, -1] = phi1(t[n+1])
        for i in range(1, nx-1):
            # [Математика] разностные аппроксимации 2-го и 1-го производных:
            # u_xx ≈ (u_{i+1} - 2u_i + u_{i-1}) / h²
            # u_x ≈ (u_{i+1} - u_{i-1}) / (2h)
            term_xx = (U[n, i+1] - 2.0*U[n, i] + U[n, i-1]) / (h*h)
            term_x  = (U[n, i+1] - U[n, i-1]) / (2.0*h)
            U[n+1, i] = 2.0*U[n, i] - U[n-1, i] + tau*tau * (2.0 * term_xx + 4.0 * term_x)
    return x, t, U

# --- Неявная конечно-разностная схема ---
# [Математика] Решается система линейных уравнений A * u^{n+1} = rhs
def implicit_scheme_from_grids(x, t, h, tau):
    nx = len(x); nt = len(t)
    U = np.zeros((nt, nx))

    # начальные условия
    for i in range(nx):
        U[0, i] = psi0(x[i])
        U[1, i] = psi0(x[i]) + tau * psi1(x[i])

    # [Математика] коэффициенты конечно-разностной аппроксимации
    alpha = tau * tau / (h * h)
    beta = tau * tau / h

    coef_center = 1.0 + 4.0 * alpha
    coef_left   = -2.0 * alpha + 2.0 * beta
    coef_right  = -2.0 * alpha - 2.0 * beta

    # [Математика] матрица трёхдиагональной СЛАУ
    m = nx - 2
    A = np.zeros((m, m))
    for i in range(m):
        A[i, i] = coef_center
        if i > 0:
            A[i, i-1] = coef_left
        if i < m-1:
            A[i, i+1] = coef_right

    # решение по слоям
    for n in range(1, nt-1):
        rhs = 2.0 * U[n, 1:-1] - U[n-1, 1:-1]
        rhs[0]  -= coef_left  * phi0(t[n+1])
        rhs[-1] -= coef_right * phi1(t[n+1])

        # [Математика] решаем A * u^{n+1} = rhs методом Гаусса
        u_interior_next = np.linalg.solve(A, rhs)

        U[n+1, 0] = phi0(t[n+1])
        U[n+1, -1] = phi1(t[n+1])
        U[n+1, 1:-1] = u_interior_next

    return x, t, U

# --- Строим аналитическое и численные решения ---
x, t, U_analytic = build_analytical_from_grids(x, t)
x, t, U_explicit = explicit_scheme_from_grids(x, t, h, tau)
x, t, U_implicit = implicit_scheme_from_grids(x, t, h, tau)

# --- Сравнение схем ---
solutions = {
    "explicit": U_explicit,
    "implicit": U_implicit
}

print("explicit: max err =", max_abs_error(U_explicit, U_analytic),
      " mean err =", mean_abs_error(U_explicit, U_analytic))
print("implicit: max err =", max_abs_error(U_implicit, U_analytic),
      " mean err =", mean_abs_error(U_implicit, U_analytic))

# --- Визуализация временной эволюции ---
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def plot_scheme_evolution_separate_windows(solutions, x, t, time_points=None):
    """
    Рисует u(x,t) для выбранных моментов времени отдельно для каждой схемы.
    """
    if time_points is None:
        time_points = [0.1, 0.2, 0.4, 0.6, 1.0]
    colors = plt.cm.viridis(np.linspace(0, 1, len(time_points)))
    schemes = list(solutions.keys())

    for scheme_name in schemes:
        plt.figure(figsize=(10, 6))
        sol = np.asarray(solutions[scheme_name])
        for i, time_point in enumerate(time_points):
            cur_t_id = int(np.abs(t - time_point).argmin())
            actual_time = t[cur_t_id]
            plt.plot(x, sol[cur_t_id, :],
                     color=colors[i],
                     label=f't={actual_time:.3f}',
                     linewidth=2)
        plt.xlabel('x', fontsize=12)
        plt.ylabel('u(x, t)', fontsize=12)
        plt.title(f'Эволюция во времени: {scheme_name}', fontsize=14)
        plt.legend(fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.xlim([x[0], x[-1]])
        plt.tight_layout()
        plt.show()

plot_scheme_evolution_separate_windows(solutions, x, t)

# --- Ошибка во времени (максимальная по x) ---
def plot_errors_from_time(solutions, analytical, t):
    plt.figure(figsize=(10,6))
    for name, sol in solutions.items():
        sol = np.asarray(sol)
        err = np.array([np.max(np.abs(sol[i,:] - analytical[i,:])) for i in range(len(t))])
        plt.plot(t, err, label=name)
    plt.xlabel('time')
    plt.ylabel('max abs error')
    plt.title('Максимальная ошибка численных схем во времени')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

plot_errors_from_time(solutions, U_analytic, t)

# --- Абсолютная ошибка на последнем временном слое ---
def plot_final_time_errors_abs(solutions, analytical, x, t):
    """
    Абсолютные ошибки |U_num - U_ana| при t = t[-1].
    """
    final_time_idx = -1
    plt.figure(figsize=(12, 8))
    for name, sol in solutions.items():
        sol = np.asarray(sol)
        abs_error = np.abs(sol[final_time_idx, :] - analytical[final_time_idx, :])
        plt.plot(x, abs_error, label=f'{name}', linewidth=2)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('Абсолютная ошибка |U_num(x) - U_ana(x)|', fontsize=12)
    plt.title(f'Абсолютные ошибки на конечном временном срезе t={t[final_time_idx]:.4f}', fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

plot_final_time_errors_abs(solutions, U_analytic, x, t)
