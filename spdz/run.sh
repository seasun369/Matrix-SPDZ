#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_tests.sh            # 运行所有 tests
#   ./run_tests.sh test1      # 只运行 test/test1.cpp 生成的可执行 test1
#   ./run_tests.sh test2 test3  # 分别运行 test2 和 test3

# 收集用户指定的目标（如果有），否则留空
declare -a TARGETS=()
if (( $# > 0 )); then
  for arg in "$@"; do
    TARGETS+=("$arg")
  done
fi

# cd build


# 遍历所有 test/*.cpp 生成的可执行，并根据 TARGETS 决定是否运行
for src in ../test/*.cpp; do
  name=$(basename "$src" .cpp)

  # 如果用户指定了 TARGETS，则只跑被指定的那些
  if (( ${#TARGETS[@]} > 0 )); then
    skip=true
    for t in "${TARGETS[@]}"; do
      if [[ "$t" == "$name" ]]; then
        skip=false
        break
      fi
    done
    $skip && continue
  fi

  echo -e "\n=== Running test: $name ==="
  ./"$name" 1 12345 &
  ./"$name" 2 12345 
done

echo -e "\nAll requested tests finished."
