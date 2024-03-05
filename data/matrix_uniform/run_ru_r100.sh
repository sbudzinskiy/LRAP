#!/bin/bash

N=(20000 10000 7500 5000 4000 3000 2500 2000 1800 1600 1400 1200 1000 900 800 700 600 500 450 400 350 300 250 200 150 125)
r=100

p=5
q=1

executable="../../build/intel-release-mkl/bin/matrix_uniform"

for n in "${N[@]}"
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

