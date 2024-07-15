# SPDX-License-Identifier: Apache-2.0

FOREACH(SVARIANT ${SVARIANT_S})
    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)

    add_executable(sqisign_test_prof_${SVARIANT} test_sqisign_prof.c)
    target_link_libraries(sqisign_test_prof_${SVARIANT} sqisign_${SVARIANT_LOWER})
    target_include_directories(sqisign_test_prof_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})

    if(ENABLE_SIGN AND ${SVARIANT_LOWER} IN_LIST SVARIANTSIGN_S)
        add_executable(sqisign_bench_${SVARIANT} bench.c)
        target_link_libraries(sqisign_bench_${SVARIANT} sqisign_${SVARIANT_LOWER})
        target_include_directories(sqisign_bench_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})
        target_compile_definitions(sqisign_bench_${SVARIANT} PUBLIC ENABLE_SIGN)

        add_executable(sqisign_test_scheme_${SVARIANT} test_sqisign.c)
        target_link_libraries(sqisign_test_scheme_${SVARIANT} sqisign_${SVARIANT_LOWER})
        target_include_directories(sqisign_test_scheme_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})
        target_compile_definitions(sqisign_test_scheme_${SVARIANT} PUBLIC ENABLE_SIGN)
    endif()

    if(NOT(${SVARIANT_LOWER} IN_LIST SVARIANTSIGN_S))
        # Executables for non-NIST-compliant variants

        add_executable(sqisign_bench_verif_kat_${SVARIANT}_smart_compressed bench_verif_kat.c)
        target_link_libraries(sqisign_bench_verif_kat_${SVARIANT}_smart_compressed sqisign_${SVARIANT_LOWER})
        target_include_directories(sqisign_bench_verif_kat_${SVARIANT}_smart_compressed PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})
        target_compile_definitions(sqisign_bench_verif_kat_${SVARIANT}_smart_compressed PRIVATE SMART_SIGNATURE)

        add_executable(sqisign_test_kat_${SVARIANT}_smart_compressed test_kat.c)
        target_link_libraries(sqisign_test_kat_${SVARIANT}_smart_compressed sqisign_${SVARIANT_LOWER}_test)
        target_include_directories(sqisign_test_kat_${SVARIANT}_smart_compressed PRIVATE ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER} ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON} ${INC_ID2ISO} ${INC_PROTOCOLS})
        target_compile_definitions(sqisign_test_kat_${SVARIANT}_smart_compressed PRIVATE SMART_SIGNATURE)
        add_test(sqisign_test_kat_${SVARIANT}_smart_compressed  sqisign_test_kat_${SVARIANT}_smart )

        add_executable(sqisign_bench_verif_kat_${SVARIANT}_smart_uncompressed bench_verif_kat.c)
        target_link_libraries(sqisign_bench_verif_kat_${SVARIANT}_smart_uncompressed sqisign_${SVARIANT_LOWER})
        target_include_directories(sqisign_bench_verif_kat_${SVARIANT}_smart_uncompressed PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})
        target_compile_definitions(sqisign_bench_verif_kat_${SVARIANT}_smart_uncompressed PRIVATE UNCOMPRESSED_SIGNATURE)

        add_executable(sqisign_test_kat_${SVARIANT}_smart_uncompressed test_kat.c)
        target_link_libraries(sqisign_test_kat_${SVARIANT}_smart_uncompressed sqisign_${SVARIANT_LOWER}_test)
        target_include_directories(sqisign_test_kat_${SVARIANT}_smart_uncompressed PRIVATE ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER} ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON} ${INC_ID2ISO} ${INC_PROTOCOLS})
        target_compile_definitions(sqisign_test_kat_${SVARIANT}_smart_uncompressed PRIVATE UNCOMPRESSED_SIGNATURE)
        add_test(sqisign_test_kat_${SVARIANT}_smart_uncompressed  sqisign_test_kat_${SVARIANT}_uncompressed )

        if(ENABLE_PARALLEL_SIGNATURE)
            add_executable(sqisign_bench_verif_kat_${SVARIANT}_parallel_uncompressed bench_verif_kat.c)
            target_link_libraries(sqisign_bench_verif_kat_${SVARIANT}_parallel_uncompressed sqisign_${SVARIANT_LOWER})
            target_include_directories(sqisign_bench_verif_kat_${SVARIANT}_parallel_uncompressed PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})
            target_compile_definitions(sqisign_bench_verif_kat_${SVARIANT}_parallel_uncompressed PRIVATE PARALLEL_SIGNATURE)

            add_executable(sqisign_test_kat_${SVARIANT}_parallel_uncompressed test_kat.c)
            target_link_libraries(sqisign_test_kat_${SVARIANT}_parallel_uncompressed sqisign_${SVARIANT_LOWER}_test)
            target_include_directories(sqisign_test_kat_${SVARIANT}_parallel_uncompressed PRIVATE ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER} ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON} ${INC_ID2ISO} ${INC_PROTOCOLS})
            target_compile_definitions(sqisign_test_kat_${SVARIANT}_parallel_uncompressed PRIVATE PARALLEL_SIGNATURE)
            add_test(sqisign_test_kat_${SVARIANT}_parallel_uncompressed  sqisign_test_kat_${SVARIANT}_parallel )

            add_executable(sqisign_bench_verif_kat_${SVARIANT}_parallel_compressed bench_verif_kat.c)
            target_link_libraries(sqisign_bench_verif_kat_${SVARIANT}_parallel_compressed sqisign_${SVARIANT_LOWER})
            target_include_directories(sqisign_bench_verif_kat_${SVARIANT}_parallel_compressed PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})
            target_compile_definitions(sqisign_bench_verif_kat_${SVARIANT}_parallel_compressed PRIVATE CPARALLEL_SIGNATURE)

            add_executable(sqisign_test_kat_${SVARIANT}_parallel_compressed test_kat.c)
            target_link_libraries(sqisign_test_kat_${SVARIANT}_parallel_compressed sqisign_${SVARIANT_LOWER}_test)
            target_include_directories(sqisign_test_kat_${SVARIANT}_parallel_compressed PRIVATE ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER} ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON} ${INC_ID2ISO} ${INC_PROTOCOLS})
            target_compile_definitions(sqisign_test_kat_${SVARIANT}_parallel_compressed PRIVATE CPARALLEL_SIGNATURE)
            add_test(sqisign_test_kat_${SVARIANT}_parallel_compressed  sqisign_test_kat_${SVARIANT}_compressed_parallel )
        endif()

    else()

        # Executables for non-NIST-compliant variants
        add_executable(sqisign_bench_verif_kat_${SVARIANT} bench_verif_kat.c)
        target_link_libraries(sqisign_bench_verif_kat_${SVARIANT} sqisign_${SVARIANT_LOWER})
        target_include_directories(sqisign_bench_verif_kat_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})

        add_executable(sqisign_test_kat_${SVARIANT} test_kat.c)
        target_link_libraries(sqisign_test_kat_${SVARIANT} sqisign_${SVARIANT_LOWER}_test)
        target_include_directories(sqisign_test_kat_${SVARIANT} PRIVATE ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER} ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON} ${INC_ID2ISO} ${INC_PROTOCOLS})

        add_test(sqisign_${SVARIANT}_KAT sqisign_test_kat_${SVARIANT})
        if(ENABLE_SIGN AND ${SVARIANT_LOWER} IN_LIST SVARIANTSIGN_S)
            add_test(sqisign_${SVARIANT}_SELFTEST sqisign_test_scheme_${SVARIANT})
        endif()
        set_tests_properties(sqisign_${SVARIANT}_KAT PROPERTIES TIMEOUT 0)

    endif()
ENDFOREACH()
