# Optimized 1D SQIsign

This is a fork of the NIST submission of SQIsign, available at https://github.com/SQISign/the-sqisign - it is advisable to read the README of this repository.

This library implements several arithmetic improvements and the 4 verification variants described in [1] for a new level-1 prime called p248, as well as an optimization of the entire protocol for the original parameter sets.

In addition to the reference (fiat-crypto) and optimized (Broadwell assembly) arithmetic options of the original library, we integrate the portable C arithmetic of Mike Scott (https://eprint.iacr.org/2024/779). We also integrate with pqm4 (https://eprint.iacr.org/2024/112), using optimised assembly code for the Cortex M4.

**NOTE:** This version corresponds to the code used to produce the results in [1].
The branch [`ePrint`](https://github.com/Crypto-TII/the-sqisign-1d/tree/ePrint)
contains code corresponding to the updated ePrint version
(https://eprint.iacr.org/2024/1563) of [1].

## Requirements

- CMake (version 3.5 or later)
- C99-compatible compiler
- Valgrind (for dynamic testing)
- Clang static analyzer (version 10 or later, for static analysis)
- GMP (version 6.1.2 or later, for signing and keygen)
- OpenMP (for parallel variants of p248)

## Build

- `mkdir -p build`
- `cd build`
- `cmake -DSQISIGN_BUILD_TYPE=<ref/broadwell/mike> -DCMAKE_BUILD_TYPE=Release ..`
- `make`

## Build options

CMake build options can be specified with `-D<BUILD_OPTION>=<VALUE>`.

In addition to the original build options, we introduce the following:

### ENABLE_SIGN

If set to `OFF`, builds only the verification of the protocol without the GMP dependency. The default value is `ON`. Note that the new prime p248 is always built with verification only, regardless of this flag.

### ENABLE_PARALLEL_SIGNATURE

If set to `ON`, cmake wil attempt to find an OpenMP installation and build the parallel verification variants for p248. Deault is `ON`


### SQISIGN_BUILD_TYPE

Specifies the build type for which SQIsign is built. The currently supported flags are:
- `ref`, which builds the plain C reference implementation.
- `broadwell`, which builds an additional implementation with GF assembly optimized code for the Intel Broadwell architecture.
- `mike`, which uses portable C arithmetic from Mike Scott's generator


## Testing

In addition to the test of the original library, tests for the new p248 verification variants can be performed by running `test/sqisign_test_kat_lvl1_p248_<variant>`, where `variant` is either `smart_compressed`, `smart_uncompressed`, `parallel_compressed` or `parallel_uncompressed`.

## Benchmarking

We have added new executables to benchmark the verification only by using the known answer test vectors. These can be performed by running `test/sqisign_bench_verif_kat_<level> <repetitions>` for the original parameter sets, and `test/sqisign_bench_verif_kat_lvl1_p248_<variant> <repetitions>` for the new parameter set with `variant` as above.


## Examples

Example instructions for reproducing our benchmarks:
```
mkdir build; cd build
cmake -DSQISIGN_BUILD_TYPE=broadwell -DCMAKE_BUILD_TYPE=Release ..
make
cd test
./sqisign_bench_verif_kat_lvl1_p248_smart_compressed 100
./sqisign_bench_verif_kat_lvl1_p248_smart_uncompressed 100
./sqisign_bench_verif_kat_lvl1_p248_parallel_compressed 100
./sqisign_bench_verif_kat_lvl1_p248_parallel_uncompressed 100
```

## pqm4

Several modifications are needed to run with pqm4, including flattening the entire directory structure and selecting only the needed source files.

This process is automated by running the script `pqm4/flatten_sources_for_pqm4.sh` script, from the root folder of the repository. It creates several folders starting with `sqisign`, which are to be copied to the `mupq/crypto_sign` folder of the pqm4 repository.

A fork of pqm4 incorporating the required changes is available at https://github.com/Crypto-TII/the-sqisign-1d-pqm4, in the `sqisign` branch.


## License

SQIsign is licensed under Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).

Third party code is used in some test and common code files:

- `src/common/aes_c.c`; MIT: "Copyright (c) 2016 Thomas Pornin <pornin@bolet.org>"
- `src/common/fips202.c`: Public Domain
- `src/common/randombytes_system.c`: MIT: Copyright (c) 2017 Daan Sprenkels <hello@dsprenkels.com>
- `apps/PQCgenKAT_sign.c`, `common/randombytes_ctrdrbg.c`, `test/test_kat.c`: by NIST (Public Domain)

## References

[1] Aardal, M. A., Adj, G., Alblooshi, A., Aranha, D. F., Canales-Martínez, I. A.,
Chávez-Saab, J., Gazzoni Filho, D. L., Reijnders, K., & Rodríguez-Henríquez, F.
Optimized One-Dimensional SQIsign Verification on Intel and Cortex-M4. IACR
Transactions on Cryptographic Hardware and Embedded Systems, 2025(1), 497-522.
https://doi.org/10.46586/tches.v2025.i1.497-522
