#!/bin/bash

d=2
N=(1025 1050 1075 1100 1150 1200 1250 1300 1400 1500 2000 2500 3000 3500 4000 4500 5000)
r=1000

executable="../../build/intel-release-mkl/bin/tensor_eye"

for n in "${N[@]}"
do
    filename="d${d}_n${n}_r${r}_ttsvd"
    params="${d} ${n} ${r}"
    for seed in {1..5}
    do
        mkdir -p "./${filename}"
	echo "d ${d} n ${n} r ${r} seed ${seed}"
	$executable $seed $params > "./${filename}/${seed}.log"
    done
done

