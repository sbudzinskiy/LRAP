#!/bin/bash

n=5000
R=(3500)

p=5
q=1

executable="../../build/intel-release-mkl/bin/matrix_uniform"

for r in "${R[@]}"
do
    filename="n${n}_r${r}_rsvd_p${p}_q${q}"
    params="${n} ${r} 2 ${p} ${q}"
    for seed in {1..10}
    do
        mkdir -p "./${filename}"
	echo "n ${n} r ${r} seed ${seed}"
	$executable $seed $params > "./${filename}/${seed}.log"
    done
done

