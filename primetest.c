#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

void reduce_modulus(mpz_t s, const mpz_t p, mpz_t q_scratch, mpz_t r_scratch) {
    uint32_t pp = mpz_get_ui(p);
    mpz_fdiv_q_2exp(q_scratch, s, pp);
    mpz_fdiv_r_2exp(r_scratch, s, pp);
    mpz_add(s, q_scratch, r_scratch);
}

int lucaslehmer(const mpz_t p)
{
    // uses lucas lehmer primality test for candidate primes 2^p - 1

    mpz_t s;
    int i;
    mpz_t modulus;
    mpz_t q_scratch, r_scratch, s2;
    int pp = mpz_get_ui(p);

    mpz_inits(q_scratch, r_scratch, NULL);

    mpz_init2(s, pp * 2);
    mpz_init2(s2, pp * 2);
    mpz_set_ui(s, 4);
    mpz_init(modulus);
    mpz_setbit(modulus, pp);
    mpz_sub_ui(modulus, modulus, 1);
    for (i = 0; i < pp - 2; i++) {
        // s = s ^ 2 - 2
        mpz_mul(s, s, s);
        // mpn_sqr(s2->_mp_d, s->_mp_d, s->_mp_size);
        mpz_sub_ui(s, s, 2);
        while (mpz_cmp(s, modulus) > 0) {
            reduce_modulus(s, p, q_scratch, r_scratch);
        }
        if (mpz_cmp(s, modulus) == 0) {
            mpz_set_ui(s, 0);
        }
    }

    mpz_clears(q_scratch, r_scratch, s, s2, NULL);
    return mpz_cmp_ui(s, 0) == 0;
}


void mersennes(unsigned long start, unsigned long stop)
{
    time_t last_output;
    mpz_t p;
    int q;
    mpz_init_set_ui(p, start);
    last_output = time(NULL);
    while (1) {
        mpz_nextprime(p, p);
        if (mpz_cmp_ui(p, stop) >= 0) {
            break;
        }
        q = lucaslehmer(p);
        if (q) {
            fputs("2^", stdout);
            mpz_out_str(stdout, 10, p);
            puts("");
        }
        else if (time(NULL) > last_output + 10) {
            fputs(" (", stdout);
            mpz_out_str(stdout, 10, p);
            puts(")");
            last_output = time(NULL);
        }
    }
}

int main(int argc, char **argv) {
    uint32_t start, stop;
    if (argc != 3) {
        printf("Usage: %s start-num stop-num\n", argv[0]);
        exit(1);
    }

    start = atoi(argv[1]);
    stop = atoi(argv[2]);

    mersennes(start, stop);

    return 0;
}
