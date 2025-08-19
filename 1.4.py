import math

# Функция создания единичной матрицы размера n x n
def create_eye_matrix(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

# Функция транспонирования матрицы (поворот по диагонали)
def transpose(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

# Функция умножения матриц A и B
def matrix_multiply(A, B):
    n = len(A)
    m = len(B[0])
    p = len(B)
    result = [[0.0 for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for k in range(p):
                result[i][j] += A[i][k] * B[k][j]
    return result

# Восстановление исходной матрицы A по её собственным значениям и векторной матрице
def create_matrix(eigen_values, eigen_vectors):
    n = len(eigen_values)
    # Создаём диагональную матрицу собственных значений (Лямбда)
    Lambda = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        Lambda[i][i] = eigen_values[i]

    # Перемножаем V * Лямбда
    VL = matrix_multiply(eigen_vectors, Lambda)
    # Перемножаем результат на V^T (транспонированную матрицу векторов)
    A_reconstructed = matrix_multiply(VL, transpose(eigen_vectors))

    return A_reconstructed

# Проверка симметричности матрицы (A[i][j] == A[j][i] для всех i, j)
def check_symmetry(A):
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i][j] - A[j][i]) > 1e-10:
                return False
    return True

# Метод вращений Якоби для поиска собственных значений и собственных векторов
def jacobi_rotation(A, precision=1e-10, max_iter=100):
    n = len(A)
    eigenvectors = create_eye_matrix(n)  # Изначально матрица собственных векторов — единичная
    iteration = 0
    iterations = []

    # Итерации вращений до достижения заданной точности
    while iteration < max_iter:
        max_off_diag = 0  # Наибольший внедиагональный элемент
        p, q = 0, 0
        # Поиск наибольшего внедиагонального элемента
        for i in range(n):
            for j in range(i + 1, n):
                if abs(A[i][j]) > max_off_diag:
                    max_off_diag = abs(A[i][j])
                    p, q = i, j

        iterations.append(iteration + 1)

        # Если все внедиагональные элементы малы — выходим
        if max_off_diag < precision:
            break

        # Вычисление угла поворота phi для зануления A[p][q]
        if A[p][p] == A[q][q]:
            phi = math.pi / 4
        else:
            phi = 0.5 * math.atan(2 * A[p][q] / (A[p][p] - A[q][q]))

        c = math.cos(phi)  # cos(phi)
        s = math.sin(phi)  # sin(phi)

        # Формируем матрицу вращения R
        R = create_eye_matrix(n)
        R[p][p] = c
        R[p][q] = -s
        R[q][p] = s
        R[q][q] = c

        # Применяем вращение: A' = R^T * A * R
        A = matrix_multiply(matrix_multiply(transpose(R), A), R)
        # Обновляем матрицу собственных векторов
        eigenvectors = matrix_multiply(eigenvectors, R)

        iteration += 1

    # Собственные значения — диагональные элементы A после итераций
    eigenvalues = [A[i][i] for i in range(n)]
    return eigenvalues, eigenvectors, iterations

# === Основной блок программы ===
try:
    # Исходная симметричная матрица A
    A = [
        [-3, -1, 3],
        [-1, 8, 1],
        [3, 1, 5]
    ]

    print("Исходная матрица A:\n")
    for row in A:
        print([round(val, 4) for val in row])

    # Проверка симметричности
    if not check_symmetry(A):
        raise ValueError("Матрица не является симметричной.")

    # Запуск метода Якоби
    eigenvalues, eigenvectors, iter_history = jacobi_rotation([row[:] for row in A])

    print("\nСобственные значения:")
    print([round(val, 6) for val in eigenvalues])

    print("\nМатрица собственных векторов:")
    for row in eigenvectors:
        print([round(val, 6) for val in row])

    print(f"\nКоличество итераций: {len(iter_history)}")

    # Восстановление матрицы из найденных собственных значений и векторов
    A_reconstructed = create_matrix(eigenvalues, eigenvectors)

    print("\nПроверка восстановления матрицы A:")
    for row in A_reconstructed:
        print([round(val, 4) for val in row])

# Вывод ошибок, если матрица не симметрична
except ValueError as e:
    print(f"Ошибка: {e}")
