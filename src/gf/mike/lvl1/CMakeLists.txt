set(SOURCE_FILES_GF_${SVARIANT_UPPER}_MIKE
    fp_p1913_32.c
    fp_p1913_64.c
    ${GFX_DIR}/fp.c
    ${GFX_DIR}/fp2.c
)

add_library(${LIB_GF_${SVARIANT_UPPER}} ${SOURCE_FILES_GF_${SVARIANT_UPPER}_MIKE})
target_include_directories(${LIB_GF_${SVARIANT_UPPER}} PRIVATE common ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} include ${PROJECT_SOURCE_DIR}/include ${INC_GF})
target_compile_options(${LIB_GF_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)
