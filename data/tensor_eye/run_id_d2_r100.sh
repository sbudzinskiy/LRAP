#!/bin/bash

d=2
N=(105 110 115 120 125 130 135 140 145 150 175 200 225 250 275 300 350 400 450 500 550 600 650 700 750 800 850)
r=100

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

