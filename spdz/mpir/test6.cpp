#include <iostream>
#include <mpir.h>

int main() {
    mpz_t *num;
    num = new mpz_t[128*128*60];
    for (int i = 0; i < 128*128*60; ++i) {
        mpz_init(num[i]);   // 初始化每个元素
        mpz_set_ui(num[i], 1234567890005844689);
        
    }
    //std::cout << "Element size: " << mpz_sizeinbase(num[0], 2) << " bits" << std::endl;
    for (int i = 0; i < 128*128*60; ++i) {
        mpz_clear(num[i]); // 清理每个元素
    }
    delete[] num; // 释放动态分配的内存
    std::cout << "end" << std::endl;
    return 0;
}

