import random
import math


def read_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    matrix = [line.strip().split() for line in lines]
    return matrix


def sum_columns(matrix):
    return [sum(map(int, col)) for col in matrix]


def mul_matrix(A, B, m):
    # sum_a,sum_b
    sum_a = []
    sum_b = []
    for i in range(len(A[0])):
        col = [row[i] for row in A]
        sum_a.append(col)
    for i in range(len(B[0])):
        col = [row[i] for row in B]
        sum_b.append(col)
    A = sum_columns(sum_a)
    B = sum_columns(sum_b)
    A = A[1:]
    B = B[1:]
    print('sum of matrix A:', A)
    print('sum of matrix B:', B)

    # 矩阵乘法
    result = [0] * (m * m)
    for i in range(m):
        for j in range(m):
            for k in range(m):
                result[i * m + j] += A[i * m + k] * B[k * m + j]
    print('sum of matrix C:', result)
    return result


def generate_c_share_matrix(mat, m, n):
    matrix = []
    k = 600000 * m * m
    for i in range(m - 1):
        row = [str(i + 1)] + [str(random.randint(1, k)) for _ in range(n)]
        matrix.append(row)
    row = [str(m)]
    for i in range(n):
        sum = 0
        for j in range(m - 1):
            sum += int(matrix[j][i + 1])
        row += [str(mat[i] - sum)]
    matrix.append(row)
    return matrix


def save_matrix_to_file(matrix, filename):
    with open(filename, 'w') as file:
        for row in matrix:
            line = ' '.join(row) + '\n'
            file.write(line)


filename1 = 'a.txt'
filename2 = 'b.txt'
matrix1 = read_file(filename1)
matrix2 = read_file(filename2)
m = len(matrix1)  # nP
n = math.ceil(math.sqrt(len(matrix1[0]) - 1))  # i.e. m
sum_result = mul_matrix(matrix1, matrix2, n)
mac_share_matrix = generate_c_share_matrix(sum_result, m, n * n)
filename3 = 'c.txt'
save_matrix_to_file(mac_share_matrix, filename3)
print("Output successfully")
