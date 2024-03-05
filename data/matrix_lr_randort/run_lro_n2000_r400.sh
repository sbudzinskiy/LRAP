#!/bin/bash

n=2000
r=400
#K=(425 450 475 500 525 550 575 600 650 700 750 800 850 900 950 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000)
K=(405 410 415 420)

p=5
q=1

executable="../../build/intel-release-mkl/bin/matrix_lr_randort"

for k in "${K[@]}"
do
    filename="n${n}_r${r}_k${k}_rsvd_p${p}_q${q}"
    params="${n} ${r} ${k} 2 ${p} ${q}"
    for seed in {1..10}
    do
        mkdir -p "./${filename}"
	echo "n ${n} r ${r} k ${k} seed ${seed}"
	$executable $seed $params > "./${filename}/${seed}.log"
    done
done

