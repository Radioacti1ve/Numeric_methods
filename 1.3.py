# Функция для проверки диагонального преобладания матрицы A
def check_matrix(A):
    n = len(A)
    flag = False
    # Проход по всем строкам
    for i in range(n):
        diagonal = abs(A[i][i])
        row_sum = sum(abs(A[i][j]) for j in range(n) if j != i)  # сумма элементов строки кроме диагонали
        if diagonal <= row_sum:
            # Если по строке не выполнено преобладание, проверяем по столбцу
            flag = False
            for i in range(n):
                diagonal = abs(A[i][i])
                column_sum = sum(abs(A[j][i]) for j in range(n) if j != i)  # сумма элементов столбца кроме диагонали
                if diagonal <= column_sum:
                    return False  # Сходимость не гарантирована
                elif diagonal > column_sum:
                    flag = True
        elif diagonal > row_sum:
            flag = True  # Сходимость возможна
    return flag  # True если хотя бы одна проверка выполнена

# Метод простых итераций для решения СЛАУ
def simple_iteration_method(A: list[list[int]], b: list[int], x0: list[int],
                            tolerance=1e-10, max_iterations=1000):
    n = len(b)
    # Строим матрицу C и вектор d для итерационного процесса
    C = [[-A[i][j] / A[i][i] if i != j else 0 for j in range(n)] for i in range(n)]
    d = [b[i] / A[i][i] for i in range(n)]

    X0 = x0.copy()  # Начальное приближение
    iteration = 0
    iteration_history = []
    residual_history = []

    # Итерационный процесс
    while iteration < max_iterations:
        # Вычисляем новое приближение
        X = [sum(C[i][j] * X0[j] for j in range(n)) + d[i] for i in range(n)]
        # Вычисляем погрешность (максимальное изменение компоненты)
        residual = max(abs(X[i] - X0[i]) for i in range(n))

        # Сохраняем итерацию и погрешность
        iteration_history.append(iteration + 1)
        residual_history.append(residual)

        # Проверка достижения нужной точности
        if residual <= tolerance:
            # Дополнительная проверка A * X ≈ b
            for i in range(n):
                calc = sum(A[i][j] * X[j] for j in range(n))
                if abs(calc - b[i]) > 1e-6:
                    raise ValueError("СЛАУ не решена точно")
            # Возвращаем результат
            return C, X, iteration + 1, iteration_history, residual_history

        X0 = X.copy()  # Переход к следующей итерации
        iteration += 1

    # Если не сошлось за max_iterations итераций
    raise Exception("Метод не сошёлся за максимальное количество итераций")

# Метод Зейделя для решения СЛАУ
def gauss_seidel(A: list[list[int]], b: list[int], x0: list[int],
                 tolerance=1e-10, max_iterations=1000):
    n = len(b)
    x = x0.copy()  # Начальное приближение
    iteration = 0
    iteration_history = []
    residual_history = []

    # Итерационный процесс
    while iteration < max_iterations:
        for i in range(n):
            # Разделяем на сумму до и после текущего элемента
            sum1 = sum(A[i][j] * x[j] for j in range(i))
            sum2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            # Вычисляем новую компоненту
            x[i] = (b[i] - sum1 - sum2) / A[i][i]

        # Вычисляем невязку (разницу между A * x и b)
        residual = [b[i] - sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
        max_residual = max(abs(r) for r in residual)

        # Сохраняем итерацию и невязку
        iteration_history.append(iteration + 1)
        residual_history.append(max_residual)

        # Проверка достижения нужной точности
        if max_residual < tolerance:
            return x, iteration + 1, iteration_history, residual_history

        iteration += 1

    # Если не сошлось за max_iterations итераций
    raise Exception("Метод не сошёлся за максимальное количество итераций")

def validate_solution(A, x, b):
    computed_b = [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]
    print("Проверка перемножения A * x:")
    for i in range(len(b)):
        print(f"Строка {i + 1}: {computed_b[i]:.6f} ≈ {b[i]} {True if abs(computed_b[i] - b[i]) < 1e-6 else False}")

# Задание СЛАУ
A = [
    [-14, 6, 1, -5],
    [-6, 27, 7, -6],
    [7, -5, -23, -8],
    [3, -8, -7, 26]
]
B = [95, -41, 69, 27]
x0 = [0, 0, 0, 0]  # Начальное приближение

print("Матрица A:", A)
print("Вектор B:", B)

# Решение методом простых итераций
C, solution_s_i, iterations_s_i, s_i_iter_history, s_i_res_history = simple_iteration_method(A, B, x0)

# Проверка условия сходимости
Chek = check_matrix(A)
if(not Chek):
    print('Предупреждение: матрица не имеет диагонального преобладания. Сходимость не гарантирована')
else:
    print("Условие сходимости метода простых итераций выполнено:", Chek)

# Вывод результата для метода простых итераций
print("Метод простых итераций:")
print("Решение:", solution_s_i)
print("Количество итераций:", iterations_s_i)
validate_solution(A, solution_s_i, B)
# Решение методом Зейделя
solution_gauss_seidel, iterations_gauss_seidel, seidel_iter_history, seidel_res_history = gauss_seidel(A, B, x0)

# Вывод результата для метода Зейделя
print("\nМетод Зейделя:")
print("Решение:", solution_gauss_seidel)
print("Количество итераций:", iterations_gauss_seidel)
validate_solution(A, solution_gauss_seidel, B)
