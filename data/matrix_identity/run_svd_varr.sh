#!/bin/bash

n=1600
R=( 20 40 80 160 320 640 )

executable="../../build/intel-release-mkl/bin/matrix_identity"

for r in "${R[@]}"
do
    filename="n${n}_r${r}_svd"
    params="${n} ${r} 1"
    for seed in {1..5}
    do
        mkdir -p "./${filename}"
	echo "n ${n} seed ${seed}"
	$executable $seed $params > "./${filename}/${seed}.log"
    done
done

