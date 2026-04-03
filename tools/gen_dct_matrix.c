/* Generate Q14 orthonormal DCT-II matrices for N=16 and N=32.
 * C[k][n] = round(alpha(k) * cos(pi*(2n+1)*k / (2*N)) * 16384)
 * alpha(0) = sqrt(1/N), alpha(k>0) = sqrt(2/N) */
#include <stdio.h>
#include <math.h>

void gen_matrix(int N) {
    printf("static const int16_t DCT%d_MAT[%d][%d] = {\n", N, N, N);
    for (int k = 0; k < N; k++) {
        double alpha = (k == 0) ? sqrt(1.0 / N) : sqrt(2.0 / N);
        printf("    {");
        for (int n = 0; n < N; n++) {
            double val = alpha * cos(M_PI * (2*n + 1) * k / (2.0 * N)) * 16384.0;
            int ival = (int)round(val);
            printf("%6d", ival);
            if (n < N - 1) printf(",");
        }
        printf("},  /* k=%d */\n", k);
    }
    printf("};\n\n");
}

int main(void) {
    gen_matrix(16);
    gen_matrix(32);
    return 0;
}
