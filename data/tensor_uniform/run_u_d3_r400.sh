#!/bin/bash

d=3
#N=(405 410 415 420 425 430 435 440 445 450 475 500 525 550 575 600 625 650 675 700 750 800 850)
N=(800 850)
r=400

executable="../../build/intel-release-mkl/bin/tensor_uniform"

echo "d ${d} n 750 r ${r} seed 10"
$executable 10 3 750 400 > "./d${d}_n750_r${r}_ttsvd/10.log"

for n in "${N[@]}"
do
    filename="d${d}_n${n}_r${r}_ttsvd"
    params="${d} ${n} ${r}"
    for seed in {1..10}
    do
        mkdir -p "./${filename}"
	echo "d ${d} n ${n} r ${r} seed ${seed}"
	$executable $seed $params > "./${filename}/${seed}.log"
    done
done

