#include <stdio.h>
#include <mpir.h>

void print_array(const mp_limb_t *array, mp_size_t size) {
    for (int i = size - 1; i >= 0; i--) {
        printf("%lx ", array[i]); // 以十六进制格式打印 limb
    }
    printf("\n");
}

int main() {
    // 初始化两个多精度数字
    mp_limb_t a[] = {0x12345678, 0xabcdef01}; // 16进制表示，每个 limb 4 字节
    mp_limb_t b[] = {0x87654321, 0xfedcba98};
    mp_size_t n = 2; // 数组中 limb 的数量

    // 打印原始数字
    printf("原始数字 a: ");
    print_array(a, n);
    printf("原始数字 b: ");
    print_array(b, n);

    // 创建结果数组
    mp_limb_t result[n + 1]; // 结果数组比输入多一个 limb，以容纳可能的进位

    // 执行加法
    mpn_add_n(result, a, b, n);

    // 打印结果
    printf("相加结果: ");
    print_array(result, n + 1);

    return 0;
}
