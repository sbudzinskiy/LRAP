#!/bin/bash

d=2
N=(405 410 415 420 425 430 435 440 445 450 475 500 525 550 575 600 625 650 675 700 750 800 850)
r=400

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

