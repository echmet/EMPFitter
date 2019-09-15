cmake_minimum_required(VERSION 3.1)

project(ElmigParamsFitter)

set(EMPF_VERSION_MAJOR 0)
set(EMPF_VERSION_MINOR 4)
set(EMPF_VERSION_PATCH 0)

add_definitions("-DECHMET_DLL_BUILD -D_USE_MATH_DEFINES")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(EXTENDED_DEBUG OFF)

if (WIN32)
    if (MINGW)
        add_definitions("-DECHMET_COMPILER_MINGW")
    elseif (MSYS)
        add_definitions("-DECHMET_COMPILER_MSYS")
    elseif (MSVC)
        add_definitions("-DECHMET_COMPILER_MSVC")
        set(DEFAULT_SYMVER_LINK " ")
    else ()
        message(FATAL "Unsupported or mistedected compiler")
    endif ()
elseif (UNIX)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.2)
            message(FATAL_ERROR "GCC 4.9.2 or newer is required to build libElmigParamsFitter")
	    endif ()

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -Wall -Wextra -pedantic -fvisibility=hidden -g -fopenmp")
        add_definitions("-DECHMET_COMPILER_GCC_LIKE")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4)
            message(FATAL_ERROR "Clang 3.4 or newer is required to build libElmigParamsFitter")
        endif ()

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wpedantic -Wno-shadow-field-in-constructor -Wno-weak-vtables -Wno-padded -fvisibility=hidden -g -fopenmp")
        add_definitions("-DECHMET_COMPILER_GCC_LIKE")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    else ()
        message(FATAL "Unsupported or misdetected compiler")
    endif ()
endif ()

if (EXTENDED_DEBUG)
    add_definitions("-DECHMET_MOBDISSOCREGRESSOR_DEBUG_DUMP")
endif ()

include_directories(${INCLUDE_DIRECTORIES}
                    "${CMAKE_CURRENT_SOURCE_DIR}/include"
                    "${CMAKE_CURRENT_BINARY_DIR}")

if (NOT "${EIGEN_INCLUDE_DIR}" STREQUAL "")
    include_directories(${INCLUDE_DIRECTORIES}
                        SYSTEM "${EIGEN_INCLUDE_DIR}")
endif ()

if (NOT "${ECHMET_CORE_LIBS_DIR}" STREQUAL "")
    include_directories(${INCLUDE_DIRECTORIES}
                        SYSTEM "${ECHMET_CORE_LIBS_DIR}/include/ECHMET/CoreLibs")
    link_directories(${LINK_DIRECTORIES}
                     "${ECHMET_CORE_LIBS_DIR}/lib")
endif ()

if (NOT "${ECHMET_REGRESSOR_DIR}" STREQUAL "")
    include_directories(${INCLUDE_DIRECTORIES}
                        SYSTEM "${ECHMET_REGRESSOR_DIR}")
endif ()

set(libElmigParamsFitter_SRCS
    src/mobdissocregressor.cpp
    src/elmigparamsfitter.cpp
    src/parametersfixerimpl.cpp)

set(EPFITTER_LINK_LIBS
    ECHMETShared
    SysComp
    CAES
    IonProps)


CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/src/elmigparamsfitter_config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/elmigparamsfitter_config.h")
add_library(ElmigParamsFitter SHARED ${libElmigParamsFitter_SRCS})
target_include_directories(ElmigParamsFitter PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
set_target_properties(ElmigParamsFitter
                                 PROPERTIES VERSION ${EMPF_VERSION_MAJOR}.${EMPF_VERSION_MINOR}
                                            SOVERSION ${EMPF_VERSION_MAJOR}.${EMPF_VERSION_MINOR}
                                 LINK_FLAGS ${DEFAULT_SYMVER_LINK}
                      PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/echmetelmigparamsfitter.h")
target_link_libraries(ElmigParamsFitter
                      PRIVATE ${EPFITTER_LINK_LIBS})

if (WIN32)
    install(TARGETS ElmigParamsFitter
            RUNTIME DESTINATION lib
            ARCHIVE DESTINATION lib
            PUBLIC_HEADER DESTINATION "include/ECHMET")
else ()
    install(TARGETS ElmigParamsFitter
            LIBRARY DESTINATION lib
            PUBLIC_HEADER DESTINATION "include/ECHMET")
endif ()
if (WIN32)
    if (MINGW OR MSYS)
	    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libElmigParamsFitter.lib" DESTINATION lib RENAME ElmigParamsFitter.lib)
	    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libElmigParamsFitter.dll" DESTINATION lib RENAME ElmigParamsFitter.dll)
    endif ()
endif ()
