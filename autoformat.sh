#!/bin/sh

echo "Formatting..."
root=`dirname ${BASH_SOURCE[0]}`
echo $root
clang-format -i $root/*.cc
clang-format -i $root/include/*.hh
clang-format -i $root/src/*.cc
echo "Done"
