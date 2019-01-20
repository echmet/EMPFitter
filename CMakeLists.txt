cmake_minimum_required(VERSION 3.1)

add_definitions("-DECHMET_DLL_BUILD")

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

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -Wall -Wextra -pedantic -fvisibility=hidden -g")
        add_definitions("-DECHMET_COMPILER_GCC_LIKE")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4)
		message(FATAL_ERROR "Clang 3.4 or newer is required to build libElmigParamsFitter")
	endif ()

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -Weverything -Wno-c++98-compat -Wpedantic -fvisibility=hidden -g")
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
                    SYSTEM "${EIGEN_INCLUDE_DIR}"
                    SYSTEM "${ECHMET_CORE_LIBS_DIR}/include/ECHMET/CoreLibs"
                    SYSTEM "${ECHMET_REGRESSOR_DIR}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(libElmigParamsFitter_SRCS
    src/elmigparamsfitter.cpp
    src/parametersfixerimpl.cpp)

set(EPFITTER_LINK_LIBS
    ECHMETShared
    SysComp
    CAES
    IonProps)

link_directories("${ECHMET_CORE_LIBS_DIR}/lib")

add_library(ElmigParamsFitter SHARED ${libElmigParamsFitter_SRCS})
target_include_directories(ElmigParamsFitter PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
set_target_properties(ElmigParamsFitter
                      PROPERTIES VERSION 0.1
                                 SOVERSION 0.1
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
