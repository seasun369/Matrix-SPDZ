#include <stdio.h>
#include <mpir.h>

int main() {
    // 初始化MPIR
    //mpir_init();

    // 定义被除数和除数
    mp_limb_t u[3] = {0x12345678, 0x9ABCDEF0, 0x12345678}; // 被除数
    mp_limb_t v = 0xABCDEF01; // 除数

    // 计算余数
    mp_limb_t remainder = mpn_mod_1(u, 3, v);

    // 打印结果
    printf("Remainder: %lx\n", remainder);

    // 清理MPIR
    //mpir_clear();

    return 0;
}

