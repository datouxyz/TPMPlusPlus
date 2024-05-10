# CMakeList.txt: Tpm++ 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

# 将源代码添加到此项目的可执行文件。
#add_executable (Tpm++ "Tpm++.cpp" "Tpm++.h")

# collect all the source files in the src directory
file(GLOB SOURCES "*.cpp")

file(GLOB struct2xSrc "struct2x/src/json/*.cpp" "struct2x/src/thirdParty/json/*.cpp" "struct2x/src/protobuf/*.cpp" "SPMSRC/*.cpp")

include_directories(".")
include_directories("SPMSRC")
include_directories("eigen-3.4.0")
include_directories("eigen-3.4.0/unsupported")
include_directories("sitk/include/SimpleITK-2.1")

#add_definitions(-DEXPORTAPI)


if(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Debug")
    add_definitions(-D_DEBUG=1)
    link_directories("E:/sitkbuild64/ITK-build/lib/Debug")#
    link_directories("E:/sitkbuild64/SimpleITK-build/lib/Debug")
else()
    message(STATUS "release")
    add_definitions(-DNDEBUG=1)
    link_directories("E:/sitkbuild64/ITK-build/lib/MinSizeRel")
    link_directories("E:/sitkbuild64/SimpleITK-build/lib/MinSizeRel")
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    link_libraries(OpenMP::OpenMP_CXX)
endif()


if(MSVC)
  add_compile_options(/bigobj)
  
 #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

endif()

add_executable(Tpm++ ${SOURCES} ${struct2xSrc})

target_link_libraries(Tpm++ 
    SimpleITKCommon-2.1
    SimpleITKIO-2.1
    SimpleITKBasicFilters0-2.1
    SimpleITKBasicFilters1-2.1
    SimpleITK_ITKCommon-2.1
)
target_sources(Tpm++ PRIVATE stdafx.cpp)
target_precompile_headers(Tpm++ PRIVATE stdafx.h)
