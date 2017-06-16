#!/bin/sh

echo "Formatting..."
root=`dirname ${BASH_SOURCE[0]}`
echo $root
clang-format -i $root/*.cc $root/include/*.hh $root/src/*.cc
echo "Done"
