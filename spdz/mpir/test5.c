#include <stdio.h>
#include <mpir.h>

int main() {
    mpz_t a, b, sum,c, mod_result;
    mpz_inits(a, b, sum,c, mod_result, NULL);

    // Initialize variables
    mpz_set_str(a, "-15", 10); // Initialize a with a large number
    mpz_set_str(b, "9", 10); // Initialize b with another large number

    // Perform addition
    mpz_add(sum, a, b);
    mpz_sub(c,a,b);

    // Perform modulo
    mpz_mod(mod_result, a, b);

    // Print results
    gmp_printf("a - b = %Zd\n", c);
    gmp_printf("a + b = %Zd\n", sum);
    gmp_printf("a mod b = %Zd\n", mod_result);

    // Clear allocated memory
    mpz_clears(a, b, sum, mod_result, NULL);

    return 0;
}

