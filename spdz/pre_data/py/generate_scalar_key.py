import random


def generate_random_matrix(m, n):
    data = []
    for i in range(m):
        row_data = [str(i + 1)] + ['0'] * (n * n)
        num = random.randint(1, 100)
        for j in range(n):
            row_data[j * n + j + 1] = str(num)
        data.append(row_data)
    return data


def save_matrix_to_file(matrix, filename):
    with open(filename, 'w') as file:
        for row in matrix:
            line = ' '.join(row) + '\n'
            file.write(line)


m = int(input("Enter number of rows(i.e. nP):"))
n = int(input("Enter m:"))
matrix = generate_random_matrix(m, n)
filename = 'scalar_key.txt'
save_matrix_to_file(matrix, filename)
print("Output successfully")
