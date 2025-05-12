#!/bin/bash

# 清理任何先前的构建
rm -rf build
mkdir build
cd build

# 运行 CMake 配置
cmake ..

# 编译项目
make -j 96

