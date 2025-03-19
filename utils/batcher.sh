#!/bin/bash

x=`find . -type f | wc -l`
y=20
batch_size=$(expr $x / $y + 1)
echo $batch_size
for i in `seq 1 $y`
do
mkdir -p "folder$i"
find . -maxdepth 1 -type f | head -n $batch_size| xargs -i mv "{}" "folder$i"
done
