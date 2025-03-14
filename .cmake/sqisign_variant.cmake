set(MODULES_SIGN "INTBIG;KLPT;QUATERNION")

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/generic)
    set(LIB_${CCSD_NAME_UPPER} sqisign_${CCSD_NAME}_generic CACHE INTERNAL "LIB")
    set(INC_${CCSD_NAME_UPPER} ${CMAKE_CURRENT_SOURCE_DIR}/generic/include CACHE INTERNAL "LIB")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/generic)
    FOREACH(SVARIANT ${SVARIANT_S})
        string(TOUPPER ${SVARIANT} SVARIANT_UPPER)
        string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
        set(LIB_${CCSD_NAME_UPPER}_${SVARIANT_UPPER} ${LIB_${CCSD_NAME_UPPER}} CACHE INTERNAL "LIB")
        set(INC_${CCSD_NAME_UPPER}_${SVARIANT_UPPER} ${INC_${CCSD_NAME_UPPER}} CACHE INTERNAL "INC")
    ENDFOREACH()
else()
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/include)
        set(INC_${CCSD_NAME_UPPER} ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "LIB")
    endif()
    FOREACH(SVARIANT ${SVARIANT_S})
        string(TOUPPER ${SVARIANT} SVARIANT_UPPER)
        string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
        # Check if variant must be included
        if((NOT ${CCSD_NAME_UPPER} IN_LIST MODULES_SIGN) OR (ENABLE_SIGN AND ${SVARIANT_LOWER} IN_LIST SVARIANTSIGN_S))
            set(LIB_${CCSD_NAME_UPPER}_${SVARIANT_UPPER} sqisign_${CCSD_NAME}_${SVARIANT} CACHE INTERNAL "LIB")
            if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${SVARIANT})
                set(INC_${CCSD_NAME_UPPER}_${SVARIANT_UPPER} ${CMAKE_CURRENT_SOURCE_DIR}/${SVARIANT}/include CACHE INTERNAL "INC")
                if(ENABLE_SIGN AND ${SVARIANT_LOWER} IN_LIST SVARIANTSIGN_S)
                    add_compile_definitions(ENABLE_SIGN)
                endif()
                add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${SVARIANT})
                if(ENABLE_SIGN AND ${SVARIANT_LOWER} IN_LIST SVARIANTSIGN_S)
                    remove_definitions(-DENABLE_SIGN)
                endif()
            else()
                message(FATAL_ERROR "No matching implementation found for variant ${SVARIANT}")
            endif()
        endif()
    ENDFOREACH()
endif()
