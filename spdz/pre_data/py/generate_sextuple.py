import random
kk = 2305843009213693951

def generate_random_matrix(m, n):
    matrix = []
    for i in range(m):
        row = [str(i + 1)] + [str(random.randint(1, kk)) for _ in range(n)]
        matrix.append(row)
    return matrix


def read_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    matrix = [line.strip().split() for line in lines]
    return matrix


def sum_columns(matrix):
    return [sum(map(int, col)) % kk for col in matrix]


def mul_matrix(A, B, m, n, l):
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
    result = [0] * (m * l)
    for i in range(m):
        for j in range(l):
            for k in range(n):
                result[i * l + j] = (result[i * l + j] + A[i * n + k] * B[k * l + j]) % kk
    print('sum of matrix C:', result)
    return result


def generate_share_matrix(mat, m, n):
    matrix = []
    for i in range(m - 1):
        row = [str(i + 1)] + [str(random.randint(1, kk)) for _ in range(n)]
        matrix.append(row)
    row = [str(m)]
    for i in range(n):
        sum = 0
        for j in range(m - 1):
            sum = (sum + int(matrix[j][i + 1])) % kk
        row += [str((mat[i] - sum) % kk)]
    matrix.append(row)
    return matrix


def transpose_matrix(matrix):
    return [list(row) for row in zip(*matrix)]


def process_file(input_filename, output_filename, m, n):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            parts = line.strip().split()
            row_number = parts[0]  # 行号
            matrix_data = [int(x) for x in parts[1:]]  # 矩阵数据

            # 构造原始矩阵 (m * n)
            matrix = [matrix_data[i * n:(i + 1) * n] for i in range(m)]

            # 转置矩阵
            transposed_matrix = transpose_matrix(matrix)

            # 写入转置后的矩阵到输出文件
            outfile.write(row_number + ' ')  # 写入行号
            for row in transposed_matrix:
                outfile.write(' '.join(map(str, row)) + ' ')
            outfile.write('\n')


def generate_vector_key(m, n):
    matrix = []
    for i in range(m):
        row = [str(i + 1)] + [str(random.randint(1, kk)) for _ in range(n)]
        matrix.append(row)
    return matrix


def mac_matrix(matrix1, matrix2, m, n):
    transposed = []
    sum_scalar_key = []
    for i in range(len(matrix1[0])):
        col = [row[i] for row in matrix1]
        transposed.append(col)
    for i in range(len(matrix2[0])):
        col = [row[i] for row in matrix2]
        sum_scalar_key.append(col)
    transposed = sum_columns(transposed)
    transposed = transposed[1:]
    print('sum of matrix a:', transposed)
    sum_scalar_key = sum_columns(sum_scalar_key)
    sum_scalar_key = sum_scalar_key[1:]
    print('sum of vector key:', sum_scalar_key)
    result = matrix_vector_mult(transposed, sum_scalar_key, m, n)
    print('sum of vector mac:', result)
    return result


def generate_mac_share_matrix(mat, m, n):
    matrix = []
    for i in range(m - 1):
        row = [str(i + 1)] + [str(random.randint(1, kk)) for _ in range(n)]
        matrix.append(row)
    row = [str(m)]
    for i in range(n):
        sum = 0
        for j in range(m - 1):
            sum = (sum + int(matrix[j][i + 1])) % kk
        row += [str((mat[i] - sum) % kk)]
    matrix.append(row)
    return matrix


def matrix_vector_mult(matrix, vector, m, n):
    result = [0] * m
    for i in range(m):
        for j in range(n):
            result[i] = (result[i] + matrix[i * n + j] * vector[j]) % kk
    return result


def save_matrix_to_file(matrix, filename):
    with open(filename, 'w') as file:
        for row in matrix:
            line = ' '.join(row) + '\n'
            file.write(line)


nP = int(input("Enter number of rows(i.e. nP):"))
m = int(input("Enter m:"))
n = int(input("Enter n:"))
l = int(input("Enter l:"))
matrix_a = generate_random_matrix(nP, m * n)
filename1 = 'a.txt'
save_matrix_to_file(matrix_a, filename1)
matrix_b = generate_random_matrix(nP, n * l)
filename2 = 'b.txt'
save_matrix_to_file(matrix_b, filename2)
matrix_r = generate_random_matrix(nP, m * l)
filename3 = 'r.txt'
save_matrix_to_file(matrix_r, filename3)
A = read_file(filename1)
B = read_file(filename2)
R = read_file(filename3)
sum_result = mul_matrix(A, B, m, n, l)
matrix_c = generate_share_matrix(sum_result, nP, m * l)
filename4 = 'c.txt'
save_matrix_to_file(matrix_c, filename4)
C = read_file(filename4)
filename5 = 'a_t.txt'
filename6 = 'r_t.txt'
process_file(filename1, filename5, m, n)
process_file(filename3, filename6, m, l)
A_T = read_file(filename5)
R_T = read_file(filename6)
matrix_key = generate_vector_key(nP, max(m, n, l))
filename7 = 'vector_key.txt'
save_matrix_to_file(matrix_key, filename7)
key = read_file(filename7)

sum_vector_mac = mac_matrix(A, key, m, n)
mac_share_matrix = generate_mac_share_matrix(sum_vector_mac, nP, m)
filename8 = 'mac_a.txt'
save_matrix_to_file(mac_share_matrix, filename8)
sum_vector_mac = mac_matrix(B, key, n, l)
mac_share_matrix = generate_mac_share_matrix(sum_vector_mac, nP, n)
filename9 = 'mac_b.txt'
save_matrix_to_file(mac_share_matrix, filename9)
sum_vector_mac = mac_matrix(C, key, m, l)
mac_share_matrix = generate_mac_share_matrix(sum_vector_mac, nP, m)
filename10 = 'mac_c.txt'
save_matrix_to_file(mac_share_matrix, filename10)
sum_vector_mac = mac_matrix(R, key, m, l)
mac_share_matrix = generate_mac_share_matrix(sum_vector_mac, nP, m)
filename11 = 'mac_r.txt'
save_matrix_to_file(mac_share_matrix, filename11)
sum_vector_mac = mac_matrix(A_T, key, n, m)
mac_share_matrix = generate_mac_share_matrix(sum_vector_mac, nP, n)
filename12 = 'mac_a_t.txt'
save_matrix_to_file(mac_share_matrix, filename12)
sum_vector_mac = mac_matrix(R_T, key, l, m)
mac_share_matrix = generate_mac_share_matrix(sum_vector_mac, nP, l)
filename13 = 'mac_r_t.txt'
save_matrix_to_file(mac_share_matrix, filename13)
print("Output successfully")