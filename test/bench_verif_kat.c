// SPDX-License-Identifier: Apache-2.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <sig.h>
#include <api.h>
#include <encoded_sizes.h>

#if defined(TARGET_OS_UNIX) && (defined(TARGET_ARM) || defined(TARGET_ARM64) || defined(TARGET_OTHER))
#include <time.h>
#endif

#define MAX_MARKER_LEN    50

#define KAT_FILE_OPEN_ERROR    -1
#define KAT_DATA_ERROR         -3

#define LIST_SIZE    10000

#define MSG_LEN_MAX       3300    // greater or equal to the maximum value of mlen in kat file

#if defined(SMART_SIGNATURE)
#define SIGMSG_LEN_MAX    (MSG_LEN_MAX + SMART_SIGNATURE_LEN)
#elif defined(UNCOMPRESSED_SIGNATURE)
#define SIGMSG_LEN_MAX    (MSG_LEN_MAX + UNCOMPRESSED_SIGNATURE_LEN)
#elif defined(PARALLEL_SIGNATURE)
#define SIGMSG_LEN_MAX    (MSG_LEN_MAX + PARALLEL_SIGNATURE_LEN)
#elif defined(CPARALLEL_SIGNATURE)
#define SIGMSG_LEN_MAX    (MSG_LEN_MAX + COMPRESSED_PARALLEL_SIGNATURE_LEN)
#else
#define SIGMSG_LEN_MAX    (MSG_LEN_MAX + SIGNATURE_LEN)
#endif

#if (defined(TARGET_ARM) || defined(TARGET_ARM64) || defined(TARGET_S390X))
#define BENCH_UNITS "nsec"
#else
#define BENCH_UNITS "cycles"
#endif

#if (defined(TARGET_ARM) || defined(TARGET_ARM64) || defined(TARGET_S390X) || defined(TARGET_OTHER))
#define print_unit printf("nsec\n");
#else
#define print_unit printf("cycles\n");
#endif

static inline int64_t cpucycles(void);
static int cmpfunc(const void *a, const void *b);
static int bench_sig(int runs, int csv);
static int FindMarker(FILE *infile, const char *marker);
static int ReadHex(FILE *infile, unsigned char *A, int Length, char *str);

/* TODO: check */
int main(int argc, char *argv[])
{
    int rc = 0;

#ifndef NDEBUG
    fprintf(stderr, "\x1b[31mIt looks like SQIsign was compiled with assertions enabled.\n"
                    "This will severely impact performance measurements.\x1b[0m\n");
#endif

    if (argc < 2) {
        printf("One argument needed\n");
        rc = 1;
        goto end;
    }
    int runs = atoi(argv[1]);
    rc = bench_sig(runs, 0);
    end:
    return rc;
}

static inline int64_t cpucycles(void)
{
#if (defined(TARGET_AMD64) || defined(TARGET_X86))
    unsigned int hi, lo;

    asm volatile ("rdtsc" : "=a" (lo), "=d"(hi));
    return ((int64_t) lo) | (((int64_t) hi) << 32);
#elif (defined(TARGET_S390X))
    uint64_t tod;
    asm volatile("stckf %0\n" : "=Q" (tod) : : "cc");
    return (tod * 1000 / 4096);
#else
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec * 1e9 + time.tv_nsec);
#endif
}

static int cmpfunc(const void *a, const void *b)
{
    return (*(int64_t *)a - *(int64_t *)b);
}

#define BENCH_CODE_1(r) \
    cycles = 0; \
    for (i = 0; i < (r); ++i) { \
        cycles1 = cpucycles();

#define BENCH_CODE_2(r, name) \
        cycles2 = cpucycles(); \
        cycles_list[i] = (cycles2 - cycles1);\
        cycles += (cycles2 - cycles1); \
    } \
    qsort(cycles_list, (r), sizeof(int64_t), cmpfunc);\
    printf("  %-20s-> median: %2" PRId64 ", average: %2" PRId64 " (%d runs) ", (name), \
    cycles_list[(r) / 2], cycles / (r), (r)); \
    printf("%s\n", BENCH_UNITS);

static int bench_sig(int runs, int csv)
{
    unsigned char *pk_data, *sm_data, *pk, *sm, m[MSG_LEN_MAX];
    size_t mlen, *smlen, *smlen_data;
    int count, i;
    int64_t cycles, cycles1, cycles2, cycles_list[LIST_SIZE];;
    char fn_rsp[64];
    FILE *fp_rsp;

    printf("Benchmarking %s\n", CRYPTO_ALGNAME);

#if defined(SMART_SIGNATURE)
    printf("Using smart signatures.\n");
    snprintf(fn_rsp, 64, "../../KAT/PQCsignKAT_%d_%s_smart.rsp", CRYPTO_SECRETKEYBYTES, CRYPTO_ALGNAME);
#elif defined(UNCOMPRESSED_SIGNATURE)
    printf("Using uncompressed signatures.\n");
    snprintf(fn_rsp, 64, "../../KAT/PQCsignKAT_%d_%s_uncompressed.rsp", CRYPTO_SECRETKEYBYTES, CRYPTO_ALGNAME);
#elif defined(PARALLEL_SIGNATURE)
    printf("Using parallel-friendly signatures.\n");
    snprintf(fn_rsp, 64, "../../KAT/PQCsignKAT_%d_%s_parallel.rsp", CRYPTO_SECRETKEYBYTES, CRYPTO_ALGNAME);
#elif defined(CPARALLEL_SIGNATURE)
    printf("Using compressed parallel-friendly signatures.\n");
    snprintf(fn_rsp, 64, "../../KAT/PQCsignKAT_%d_%s_compressed_parallel.rsp", CRYPTO_SECRETKEYBYTES, CRYPTO_ALGNAME);
#else
    printf("Using regular signatures.\n");
    snprintf(fn_rsp, 64, "../../KAT/PQCsignKAT_%d_%s.rsp", CRYPTO_SECRETKEYBYTES, CRYPTO_ALGNAME);
#endif

    if ((fp_rsp = fopen(fn_rsp, "r")) == NULL) {
        printf("Couldn't open <%s> for read\n", fn_rsp);
        return KAT_FILE_OPEN_ERROR;
    }

    if (runs > LIST_SIZE)
        runs = LIST_SIZE;

    /* Buffers for public keys and signatures + messages */
    pk_data = calloc(runs * CRYPTO_PUBLICKEYBYTES, 1);
    sm_data = calloc(runs * SIGMSG_LEN_MAX, 1);
    smlen_data = calloc(runs, sizeof(size_t));

    /* Read data from kat */
    pk = pk_data;
    sm = sm_data;
    smlen = smlen_data;
    count = 0;
    do {
        if (FindMarker(fp_rsp, "count = ")) {
            if (count >= runs)
                break;
        } else {
            break;
        }

        if (!ReadHex(fp_rsp, pk, CRYPTO_PUBLICKEYBYTES, "pk = ")) {
            printf("ERROR: unable to read 'pk' from <%s>\n", fn_rsp);
            return KAT_DATA_ERROR;
        }

        if (FindMarker(fp_rsp, "smlen = ")) {
            i = fscanf(fp_rsp, "%zd", smlen);
        } else {
            printf("ERROR: unable to read 'smlen' from <%s>\n", fn_rsp);
            return KAT_DATA_ERROR;
        }

        if (!ReadHex(fp_rsp, sm, (int)(*smlen), "sm = ")) {
            printf("ERROR: unable to read 'sm' from <%s>\n", fn_rsp);
            return KAT_DATA_ERROR;
        }

        pk += CRYPTO_PUBLICKEYBYTES;
        sm += SIGMSG_LEN_MAX;
        smlen++;
        count++;
    } while(1);

    pk = pk_data;
    sm = sm_data;
    smlen = smlen_data;

    BENCH_CODE_1(count);
#if defined(SMART_SIGNATURE)
        sqisign_open_smart(m, &mlen, sm, *smlen, pk);
#elif defined(UNCOMPRESSED_SIGNATURE)
        sqisign_open_uncompressed(m, &mlen, sm, *smlen, pk);
#elif defined(PARALLEL_SIGNATURE)
        sqisign_open_parallel(m, &mlen, sm, *smlen, pk);
#elif defined(CPARALLEL_SIGNATURE)
        sqisign_open_cparallel(m, &mlen, sm, *smlen, pk);
#else
        sqisign_open(m, &mlen, sm, *smlen, pk);
#endif
        pk += CRYPTO_PUBLICKEYBYTES;
        sm += SIGMSG_LEN_MAX;
        smlen++;
    BENCH_CODE_2(count, "sqisign_verify");

    free(pk_data);
    free(sm_data);
    free(smlen_data);

    return 0;
}

//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
static int
FindMarker(FILE *infile, const char *marker) {
    char    line[MAX_MARKER_LEN];
    int     i, len;
    int curr_line;

    len = (int)strlen(marker);
    if ( len > MAX_MARKER_LEN - 1 ) {
        len = MAX_MARKER_LEN - 1;
    }

    for ( i = 0; i < len; i++ ) {
        curr_line = fgetc(infile);
        line[i] = curr_line;
        if (curr_line == EOF ) {
            return 0;
        }
    }
    line[len] = '\0';

    while ( 1 ) {
        if ( !strncmp(line, marker, len) ) {
            return 1;
        }

        for ( i = 0; i < len - 1; i++ ) {
            line[i] = line[i + 1];
        }
        curr_line = fgetc(infile);
        line[len - 1] = curr_line;
        if (curr_line == EOF ) {
            return 0;
        }
        line[len] = '\0';
    }

    // shouldn't get here
    return 0;
}

//
// ALLOW TO READ HEXADECIMAL ENTRY (KEYS, DATA, TEXT, etc.)
//
static int
ReadHex(FILE *infile, unsigned char *A, int Length, char *str) {
    int         i, ch, started;
    unsigned char   ich;

    if ( Length == 0 ) {
        A[0] = 0x00;
        return 1;
    }
    memset(A, 0x00, Length);
    started = 0;
    if ( FindMarker(infile, str) )
        while ( (ch = fgetc(infile)) != EOF ) {
            if ( !isxdigit(ch) ) {
                if ( !started ) {
                    if ( ch == '\n' ) {
                        break;
                    } else {
                        continue;
                    }
                } else {
                    break;
                }
            }
            started = 1;
            if ( (ch >= '0') && (ch <= '9') ) {
                ich = ch - '0';
            } else if ( (ch >= 'A') && (ch <= 'F') ) {
                ich = ch - 'A' + 10;
            } else if ( (ch >= 'a') && (ch <= 'f') ) {
                ich = ch - 'a' + 10;
            } else { // shouldn't ever get here
                ich = 0;
            }

            for ( i = 0; i < Length - 1; i++ ) {
                A[i] = (A[i] << 4) | (A[i + 1] >> 4);
            }
            A[Length - 1] = (A[Length - 1] << 4) | ich;
        } else {
        return 0;
    }

    return 1;
}
