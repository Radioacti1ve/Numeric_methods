import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def tridiagonal_solve(a_diag, b_sub, c_sup, d): #решает трёхдиагональную СЛАУ методом Томаса (прогонкой)
    n = len(a_diag)
    a = a_diag.astype(float).copy()
    b = b_sub.astype(float).copy()
    c = c_sup.astype(float).copy()
    d = d.astype(float).copy()

    for i in range(1, n):
        w = b[i-1] / a[i-1]
        a[i] = a[i] - w * c[i-1]
        d[i] = d[i] - w * d[i-1]
    x = np.zeros(n, dtype=float)
    x[-1] = d[-1] / a[-1]
    for i in range(n-2, -1, -1):
        x[i] = (d[i] - c[i] * x[i+1]) / a[i]
    return x

a = 1.0
x_begin = 0.0
x_end = 1.0
t_begin = 0.0
t_end = 1.0
h = 0.01
sigma = 0.45  # sigma = tau * a / h^2 => tau = sigma * h^2 / a, где tau - шаг по времени
# для явной схемы стабильность требует sigma <= 0.5
def f_source(x):
    return np.sin(np.pi * x)

def phi_0(t):
    return 0.0

def phi_1(t):
    return 0.0

def psi(x):
    return 0.0

def analytical_solution_fn(x, t):
    return (1.0 / (math.pi**2)) * (1.0 - math.exp(- (math.pi**2) * t)) * math.sin(math.pi * x)

def get_analytical_solution(x_range, t_range, h, sigma, a=a):
    tau = sigma * h**2 / a
    x = np.arange(x_range[0], x_range[1] + h/2, h)
    t = np.arange(t_range[0], t_range[1] + tau/2, tau)
    U = np.zeros((len(t), len(x)))
    for it in range(len(t)):
        for ix in range(len(x)):
            U[it, ix] = analytical_solution_fn(x[ix], t[it])
    return x, t, U

def max_abs_error(A, B):
    assert A.shape == B.shape
    return np.abs(A - B).max()

def mean_abs_error(A, B):
    assert A.shape == B.shape
    return np.abs(A - B).mean()

# реализует метод конечных разностей (явный)
def explicit_finite_difference_method(x_range, t_range, h, sigma, a=a, phi_0=phi_0, phi_1=phi_1, psi=psi):
    tau = sigma * h**2 / a # шаг по времени
    x = np.arange(x_range[0], x_range[1] + h/2, h) # пространственная сетка
    t = np.arange(t_range[0], t_range[1] + tau/2, tau) # временная сетка

    Nx = len(x)
    Nt = len(t)
    res = np.zeros((Nt, Nx))
    for i in range(Nx): # обрабатываем начальное условие u(x, 0) = 0
        res[0, i] = psi(x[i])
    f = f_source(x)

    for n in range(1, Nt):
        res[n, 0] = phi_0(t[n]) # u(0, t) = 0
        res[n, -1] = phi_1(t[n]) # u(1, t) = 0
        for i in range(1, Nx - 1):
            res[n, i] = (
                sigma * res[n-1, i-1]
                + (1 - 2 * sigma) * res[n-1, i]
                + sigma * res[n-1, i+1]
                + tau * f[i]
            )
    return x, t, res

# реализует метод конечных разностей (неявный)
def implicit_finite_difference_method(x_range, t_range, h, sigma, a=a, phi_0=phi_0, phi_1=phi_1, psi=psi):
    tau = sigma * h**2 / a
    x = np.arange(x_range[0], x_range[1] + h/2, h)
    t = np.arange(t_range[0], t_range[1] + tau/2, tau)

    Nx = len(x)
    Nt = len(t)
    res = np.zeros((Nt, Nx))
    for i in range(Nx):
        res[0, i] = psi(x[i])
    f = f_source(x)

    n_in = Nx - 2
    if n_in > 0:
        a_diag = np.full(n_in, 1 + 2 * sigma)
        b_sub = np.full(n_in - 1, -sigma)
        c_sup = np.full(n_in - 1, -sigma)

    for n in range(1, Nt):
        # правая часть
        b = res[n-1, 1:-1].copy() + tau * f[1:-1]
        b[0] += sigma * phi_0(t[n])
        b[-1] += sigma * phi_1(t[n])

        res[n, 0] = phi_0(t[n]) # u(0, t) = 0
        res[n, -1] = phi_1(t[n]) # u(1, t) = 0
        if n_in > 0:
            res[n, 1:-1] = tridiagonal_solve(a_diag, b_sub, c_sup, b)
    return x, t, res

# метод Кранка-Никольсона
def crank_nicolson_method(x_range, t_range, h, sigma, a=a, phi_0=phi_0, phi_1=phi_1, psi=psi, theta=0.5):
    tau = sigma * h**2 / a
    x = np.arange(x_range[0], x_range[1] + h/2, h)
    t = np.arange(t_range[0], t_range[1] + tau/2, tau)

    Nx = len(x)
    Nt = len(t)
    res = np.zeros((Nt, Nx))
    for i in range(Nx):
        res[0, i] = psi(x[i])
    f = f_source(x)

    n_in = Nx - 2
    if n_in > 0:
        a_diag = np.full(n_in, 1 + 2 * sigma * theta)
        b_sub = np.full(n_in - 1, -sigma * theta)
        c_sup = np.full(n_in - 1, -sigma * theta)

    for n in range(1, Nt):
        lap_prev = np.array([res[n-1, i-1] - 2*res[n-1, i] + res[n-1, i+1] for i in range(1, Nx-1)])
        # u_{j-1}^k - 2u_j^k + u_{j+1}^k
        b = (res[n-1, 1:-1].copy() # u_j^k
              + sigma * (1 - theta) * lap_prev # + 0.5σ(u_{j-1}^k - 2u_j^k + u_{j+1}^k)
              + tau * ( (1-theta) * f[1:-1] + theta * f[1:-1] )) # + τf_j
        b[0] += sigma * theta * phi_0(t[n]) # u0^{k+1} - известна, поэтому переносим в правую часть
        b[-1] += sigma * theta * phi_1(t[n]) # последняя внутр точка правой границы тоже известна
        res[n, 0] = phi_0(t[n])
        res[n, -1] = phi_1(t[n])
        if n_in > 0:
            res[n, 1:-1] = tridiagonal_solve(a_diag, b_sub, c_sup, b)
    return x, t, res

x_ana, t_ana, analytical = get_analytical_solution(
    x_range=(x_begin, x_end),
    t_range=(t_begin, t_end),
    h=h,
    sigma=sigma,
    a=a,
)

solutions = {}
x_exp, t_exp, explicit_sol = explicit_finite_difference_method(
    x_range=(x_begin, x_end),
    t_range=(t_begin, t_end),
    h=h,
    sigma=sigma,
)
solutions["Явное"] = explicit_sol

x_imp, t_imp, implicit_sol = implicit_finite_difference_method(
    x_range=(x_begin, x_end),
    t_range=(t_begin, t_end),
    h=h,
    sigma=sigma,
)
solutions["Неявное"] = implicit_sol

x_cn, t_cn, cn_sol = crank_nicolson_method(
    x_range=(x_begin, x_end),
    t_range=(t_begin, t_end),
    h=h,
    sigma=sigma,
    theta=0.5
)
solutions["Кранк-Николсон"] = cn_sol

print("shapes:", explicit_sol.shape, analytical.shape)
print("Явное: max, mean error =",
      max_abs_error(explicit_sol, analytical),
      mean_abs_error(explicit_sol, analytical))
print("Неявное: max, mean error =",
      max_abs_error(implicit_sol, analytical),
      mean_abs_error(implicit_sol, analytical))
print("Кранк-Николсон: max, mean error =",
      max_abs_error(cn_sol, analytical),
      mean_abs_error(cn_sol, analytical))

def plot_results(solutions, time, x, t):
    cur_t_id = int(np.abs(t - time).argmin())
    plt.figure(figsize=(10,6))
    for name, sol in solutions.items():
        plt.plot(x, sol[cur_t_id], label=name)
    plt.plot(x, analytical[cur_t_id], '--', label='Аналит. решение')
    plt.xlabel('x')
    plt.ylabel(f'u(x, t={t[cur_t_id]:.4f})')
    plt.title(f'Сравнение численных и аналитического решения при t={t[cur_t_id]:.4f}')
    plt.legend()
    plt.grid()
    plt.show()

def plot_errors_from_time(solutions, analytical, t):
    plt.figure(figsize=(10,6))
    for name, sol in solutions.items():
        err = np.array([max_abs_error(sol[i], analytical[i]) for i in range(len(t))])
        plt.plot(t, err, label=name)
    plt.xlabel('time')
    plt.ylabel('max abs error')
    plt.title('Максимальная ошибка численных схем во времени')
    plt.legend()
    plt.grid()
    plt.show()


plot_results(solutions, time=0.2, x=x_ana, t=t_ana)
plot_errors_from_time(solutions, analytical, t_ana)