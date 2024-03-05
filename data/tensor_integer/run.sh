#!/bin/bash

d=3
n=100
r=5

method_lr=1

executable="../../build/intel-release-mklseq/bin/tensor_integer"
filename="d${d}_n${n}_r${r}"
params="${d} ${n} ${r}"
if [ $method_lr = 1 ]; then
    filename="${filename}_ttsvd"
fi

for seed in {1..20}
do
    mkdir -p "./${filename}"
    echo "seed ${seed}"
    $executable $seed $params > "./${filename}/${seed}.log"
done
