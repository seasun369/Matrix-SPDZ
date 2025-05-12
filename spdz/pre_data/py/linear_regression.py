import random


def generate_linear_data_with_noise(m, n, np=2):
    """
    生成线性数据，并将样本划分为np个参与方。
    每个参与方的数据是每列x_i和y的部分数据之和。

    参数：
    m：样本数量（行数）。
    n：每行数据的特征数量（不包括编号列）。
    np：参与方数量，默认为2。
    """
    # 随机生成系数w_0, w_1, ..., w_n，范围为1到50000
    weights = [random.randint(1, 50001) for _ in range(n)]  # 系数w是从1到50000的随机整数
    # print(weights)

    # 生成m个样本，每个样本有n-1个输入特征，特征值随机选择1到50000之间的整数
    X_matrix = [[random.randint(1, 50001) for _ in range(n - 1)] for _ in range(m)]  # X的特征矩阵

    # 计算 y = w_n * x_n + ... + w_1 * x_1 + w_0 + 噪声（高斯白噪声）
    y_values = custom_dot(X_matrix, weights, m, n)

    X_matrix1 = [[] for _ in range(m)]
    X_matrix2 = [[] for _ in range(m)]
    y_values1 = []
    y_values2 = []

    # 分配样本到参与方
    kk = 2305843009213693951
    for i in range(m):
        # 生成参与方 1 和参与方 2 的x_i值
        for j in range(n - 1):
            x_1 = random.randint(1, kk)  # 参与方1的份额
            x_2 = (X_matrix[i][j] - x_1) % kk  # 参与方2的份额
            X_matrix1[i].append(x_1)  # 添加到对应行
            X_matrix2[i].append(x_2)  # 添加到对应行

        # 分配 y 值
        y_1 = random.randint(1, kk)
        y_2 = (y_values[i] - y_1) % kk
        y_values1.append(y_1)
        y_values2.append(y_2)


    # 格式化结果为输出矩阵
    output_matrix = []

    # 添加参与方 1 的数据
    row=[str(1)]
    for i in range(m):
        row = row + [str(X_matrix1[i][j]) for j in range(n - 1)] + [str(y_values1[i])]
    output_matrix.append(row)

    # 添加参与方 2 的数据
    row = [str(2)]
    for i in range(m):
        row = row + [str(X_matrix2[i][j]) for j in range(n - 1)] + [str(y_values2[i])]
    output_matrix.append(row)

    return output_matrix

def custom_dot(X_matrix, weights, m, n):
    result = []
    for i in range(m):
        sum_value = weights[0]  # 初始化为w_0
        for j in range(n - 1):
            sum_value += X_matrix[i][j] * weights[j + 1]
        # 添加高斯白噪声
        noise = random.gauss(0, 1000)  # 生成高斯白噪声（均值0，标准差1000）
        noise = round(noise)  # 将噪声四舍五入为整数
        # print(noise)
        sum_value += noise
        result.append(sum_value)
    return result

def save_matrix_to_txt(matrix, filename):
    """
    将生成的数据保存到文本文件，每行的数据用空格分隔。
    """
    with open(filename, 'w') as file:
        for row in matrix:
            line = ' '.join(row) + '\n'
            file.write(line)


# 输入样本数量、特征数量和参与方数量
m = int(input("Enter the number of samples m (rows): "))  # 样本数量
n = int(input("Enter the number of features n (columns,w_0,...,w_n): "))  # 特征列数（x列数加上y列）
np = 2  # 参与方数量，默认为2

# 生成线性数据
matrix = generate_linear_data_with_noise(m, n + 1, np)  # n+1表示包含y的列数

# 保存到txt文件, x_1,...,x_n,y的顺序保存m个样本
txt_filename = 'training_data_noise.txt'
save_matrix_to_txt(matrix, txt_filename)

print(f"Data saved to {txt_filename}")
