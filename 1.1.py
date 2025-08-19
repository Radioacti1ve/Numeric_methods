def determinant(A):
    n = len(A)
    det = 1
    for i in range(n):
        det *= A[i][i]
    return det

def matrix_multiply(A, B):
    n = len(A)
    m = len(B[0])
    result = [[0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            result[i][j] = sum(A[i][k] * B[k][j] for k in range(len(B)))
    return result

def lu_decomposition(A):
    n = len(A)
    for k in range(n):
        for i in range(k + 1, n):
            A[i][k] /= A[k][k]
            for j in range(k + 1, n):
                A[i][j] -= A[i][k] * A[k][j]
    return A

def solve_lu(A, b):
    n = len(A)
    
    y = [0] * n
    for i in range(n):
        y[i] = b[i] - sum(A[i][j] * y[j] for j in range(i))
    
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]
    
    return x

def inverse_matrix(A):
    n = len(A)
    inv = [[0] * n for _ in range(n)]
    for i in range(n):
        e = [0] * n
        e[i] = 1
        inv_col = solve_lu(A, e)
        for j in range(n):
            inv[j][i] = inv_col[j]
    return inv

def check_solution(A, x, b):
    n = len(A)
    Ax = [sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
    return all(abs(Ax[i] - b[i]) < 1e-6 for i in range(n))

A = [
    [-9, 8, 8, 6],
    [-7, -9, 5, 4],
    [-3, -1, 8, 0],
    [3, -1, -4, -5]
]
b = [-81, -50, -69, 48]

A_lu = lu_decomposition([row[:] for row in A])
x = solve_lu(A_lu, b)
det_A = determinant(A_lu)
inv_A = inverse_matrix(A_lu)
is_correct = check_solution(A, x, b)

print("Решение СЛАУ (x):", x)
print("\nОпределитель матрицы A:", det_A)
print("\nОбратная матрица A:")
for row in inv_A:
    print(row)
print("\nПроверка решения (Ax = b):", is_correct)
