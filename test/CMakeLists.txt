cmake_minimum_required(VERSION 3.6)

set(CMAKE_CONFIGURATION_TYPES "Debug" "Release")

set(PROJECTNAME TestRangeCoder)
project(${PROJECTNAME})

set(HEADERS "../cpprcoder.h;slz4.h;../cppans.h;../cppase.h;../blksort.h")
set(SOURCES "main.cpp;slz4.cpp;../blksort.cpp")

include_directories(AFTER ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(AFTER ${CMAKE_CURRENT_SOURCE_DIR}/zstd/include)
link_directories("zlib")
link_directories("lz4")
link_directories("zstd/static")

source_group("include" FILES ${HEADERS})
source_group("src" FILES ${SOURCES})

set(FILES ${HEADERS} ${SOURCES})

set(OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${OUTPUT_DIRECTORY}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${OUTPUT_DIRECTORY}")

add_library(ZSTD MODULE IMPORTED)
set_target_properties(ZSTD PROPERTIES IMPORTED_LOCATION_DEBUG "libzstd_staticd.lib" IMPORTED_LOCATION_RELEASE "libzstd_static.lib")

add_executable(${PROJECTNAME} ${FILES})

set(CXX_FLAGS_SIMD "")
if(MSVC)
    if(USE_SIMD)
        set(CXX_FLAGS_SIMD "/DCPPRCODER_USE_SIMD")
    endif()

    set(DEFAULT_CXX_FLAGS "/DWIN32 /D_WINDOWS /D_MBCS /W4 /WX- /nologo /fp:precise /std:c++20 /arch:AVX /Zc:wchar_t /TP /Gd ${CXX_FLAGS_SIMD}")
    if("1800" VERSION_LESS MSVC_VERSION)
        set(DEFAULT_CXX_FLAGS "${DEFAULT_CXX_FLAGS} /EHsc")
    endif()

    set(CMAKE_CXX_FLAGS "${DEFAULT_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "/D_DEBUG /MTd /Zi /Ob0 /Od /RTC1 /Gy /GR- /GS /Gm-")
    set(CMAKE_CXX_FLAGS_RELEASE "/MT /O2 /GL /GR- /DNDEBUG")
    target_link_libraries(${PROJECTNAME} "zlibstat.lib")
    target_link_libraries(${PROJECTNAME} "liblz4_static.lib")
    target_link_libraries(${PROJECTNAME} ZSTD)

elseif(UNIX)
    if(USE_SIMD)
        set(CXX_FLAGS_SIMD "-DCPPRCODER_USE_SIMD")
    endif()
    set(DEFAULT_CXX_FLAGS "-Wall -O2 -std=c++20 -march=x86-64-v3 ${CXX_FLAGS_SIMD}")
    set(CMAKE_CXX_FLAGS "${DEFAULT_CXX_FLAGS}")
    target_link_libraries(${PROJECTNAME} "z")
elseif(APPLE)
endif()

set_target_properties(${PROJECTNAME} PROPERTIES OUTPUT_NAME_DEBUG "${PROJECTNAME}" OUTPUT_NAME_RELEASE "${PROJECTNAME}")
