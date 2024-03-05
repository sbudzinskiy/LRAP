#!/bin/bash

n=2000
r=800
#B=(1 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000)
B=(10 20 30 40 50 60 70 80 90)

p=5
q=1

executable="../../build/intel-release-mkl/bin/matrix_banded_uniform"

for b in "${B[@]}"
do
    filename="n${n}_r${r}_b${b}_rsvd_p${p}_q${q}"
    params="${n} ${r} ${b} 2 ${p} ${q}"
    for seed in {1..10}
    do
        mkdir -p "./${filename}"
	echo "n ${n} r ${r} b ${b} seed ${seed}"
	$executable $seed $params > "./${filename}/${seed}.log"
    done
done

