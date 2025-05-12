#include <stdio.h>
#include <mpir.h>

int main() {
    mpz_t n;
    mpz_init(n);
    unsigned long int a;

    // 设置大整数 n 的值
    mpz_set_str(n, "123456789012345", 10);

    // 将大整数 n 转换为 unsigned long int 类型
    a = mpz_get_ui(n);

    // 打印 unsigned long int 类型的值
    printf("%lu\n", a); // 使用 %lu 表示 unsigned long int 类型

    mpz_clear(n);
    return 0;
}

