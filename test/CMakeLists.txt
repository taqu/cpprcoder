cmake_minimum_required(VERSION 3.6)

set(CMAKE_CONFIGURATION_TYPES "Debug" "Release")

set(ProjectName TestRangeCoder)
project(${ProjectName})

set(SOURCES "main.cpp")

include_directories(AFTER ${CMAKE_CURRENT_SOURCE_DIR})
link_directories("zlib")
link_directories("lz4")

source_group("include" FILES ${HEADERS})
source_group("src" FILES ${SOURCES})

set(FILES ${HEADERS} ${SOURCES})

set(OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${OUTPUT_DIRECTORY}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${OUTPUT_DIRECTORY}")

add_executable(${ProjectName} ${FILES})

set(CXX_FLAGS_SIMD "")
if(MSVC)
    if(USE_SIMD)
        set(CXX_FLAGS_SIMD "/arch:AVX")
    endif()

    set(DEFAULT_CXX_FLAGS "/DWIN32 /D_WINDOWS /D_MBCS /W4 /WX- /nologo /fp:precise /Zc:wchar_t /TP /Gd ${CXX_FLAGS_SIMD}")
    if("1800" VERSION_LESS MSVC_VERSION)
        set(DEFAULT_CXX_FLAGS "${DEFAULT_CXX_FLAGS} /EHsc")
    endif()

    set(CMAKE_CXX_FLAGS "${DEFAULT_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "/D_DEBUG /DMU3_EXPORTS /MTd /Zi /Ob0 /Od /RTC1 /Gy /GR- /GS /Gm-")
    set(CMAKE_CXX_FLAGS_RELEASE "/DMU3_EXPORTS /MT /O2 /GL /GR- /DNDEBUG")
    target_link_libraries(${ProjectName} "zlibstat.lib")
    target_link_libraries(${ProjectName} "liblz4_static.lib")

elseif(UNIX)
    if(USE_SIMD)
        set(CXX_FLAGS_SIMD "-march=native")
    endif()
    set(DEFAULT_CXX_FLAGS "-Wall -O2 -std=c++11 ${CXX_FLAGS_SIMD}")
    set(CMAKE_CXX_FLAGS "${DEFAULT_CXX_FLAGS}")
    target_link_libraries(${ProjectName} "z")
elseif(APPLE)
endif()

set_target_properties(${ProjectName} PROPERTIES OUTPUT_NAME_DEBUG "${ProjectName}" OUTPUT_NAME_RELEASE "${ProjectName}")
