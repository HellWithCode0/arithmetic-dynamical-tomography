/*
 * cycles.c -- periodic structure of f_c(x) = x^2 + c on Z/nZ.
 *
 * Usage:
 *   cycles primes PMIN PMAX c1 [c2 ...]   one line per (prime p, probe c)
 *   cycles moduli NMIN NMAX c1 [c2 ...]   same for every modulus n in range
 *   cycles list n1,n2,... c1 [c2 ...]      explicit list of moduli
 *
 * Output line:  n c N_per | len:count len:count ...
 * where N_per is the number of periodic points and the list is the
 * multiset of cycle lengths (the isomorphism class of the periodic Z-set).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int is_prime(uint64_t n) {
    if (n < 2) return 0;
    for (uint64_t d = 2; d * d <= n; d++) if (n % d == 0) return 0;
    return 1;
}

static int cmp_u32(const void *a, const void *b) {
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}

static uint32_t *state, *pos, *path, *lens;

static void run(uint64_t n, long long c) {
    long long cm = ((c % (long long)n) + (long long)n) % (long long)n;
    memset(state, 0, n * sizeof(uint32_t));
    size_t nl = 0;
    uint64_t nper = 0;
    for (uint64_t s = 0; s < n; s++) {
        if (state[s]) continue;
        uint64_t x = s; uint32_t len = 0;
        while (!state[x]) {
            state[x] = 1; pos[x] = len; path[len++] = (uint32_t)x;
            x = (n < (1ULL << 31)) ? (x * x + (uint64_t)cm) % n
                                   : (uint64_t)(((unsigned __int128)x * x + (uint64_t)cm) % n);
        }
        if (state[x] == 1) {           /* closed a new cycle */
            uint32_t L = len - pos[x];
            lens[nl++] = L; nper += L;
        }
        for (uint32_t i = 0; i < len; i++) state[path[i]] = 2;
    }
    qsort(lens, nl, sizeof(uint32_t), cmp_u32);
    printf("%llu %lld %llu |", (unsigned long long)n, c, (unsigned long long)nper);
    for (size_t i = 0; i < nl;) {
        size_t j = i; while (j < nl && lens[j] == lens[i]) j++;
        printf(" %u:%zu", lens[i], j - i);
        i = j;
    }
    printf("\n");
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: see source\n"); return 1; }
    int mode_primes = !strcmp(argv[1], "primes");
    int mode_list = !strcmp(argv[1], "list");
    uint64_t lo = 0, hi = 0; int ci = 4;
    uint64_t *list = NULL; size_t nlist = 0;
    if (mode_list) {
        char *s = strdup(argv[2]); list = malloc(sizeof(uint64_t) * 100000);
        for (char *t = strtok(s, ","); t; t = strtok(NULL, ",")) list[nlist++] = strtoull(t, 0, 10);
        for (size_t i = 0; i < nlist; i++) if (list[i] > hi) hi = list[i];
        ci = 3;
    } else { lo = strtoull(argv[2], 0, 10); hi = strtoull(argv[3], 0, 10); }
    state = malloc((hi + 1) * sizeof(uint32_t)); pos = malloc((hi + 1) * sizeof(uint32_t));
    path = malloc((hi + 1) * sizeof(uint32_t)); lens = malloc((hi + 1) * sizeof(uint32_t));
    if (mode_list) {
        for (size_t i = 0; i < nlist; i++) for (int k = ci; k < argc; k++) run(list[i], atoll(argv[k]));
    } else {
        for (uint64_t n = lo < 2 ? 2 : lo; n <= hi; n++) {
            if (mode_primes && !is_prime(n)) continue;
            for (int k = ci; k < argc; k++) run(n, atoll(argv[k]));
        }
    }
    return 0;
}
