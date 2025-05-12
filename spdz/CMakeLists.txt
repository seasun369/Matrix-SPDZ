cmake_minimum_required(VERSION 3.14)
project(MulTest)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)

# 开启 OpenMP 支持
find_package(OpenMP REQUIRED)
# 查找并链接 MPIR 库（假设安装了 MPIR 库）·
find_library(MPIR_LIB mpir)
# 收集所有 test/*.cpp
file(GLOB TEST_SOURCES
    "${CMAKE_SOURCE_DIR}/test/*.cpp"
)

# 为每个源文件单独生成一个可执行
foreach(src IN LISTS TEST_SOURCES)
    # NAME_WE: 去掉路径与扩展名后的文件名
    get_filename_component(test_name ${src} NAME_WE)

    add_executable(${test_name} ${src})

    # 链接 OpenMP 和 pthread
    target_link_libraries(${test_name} PUBLIC
        OpenMP::OpenMP_CXX
        pthread
        ${MPIR_LIB}
    )
endforeach()
