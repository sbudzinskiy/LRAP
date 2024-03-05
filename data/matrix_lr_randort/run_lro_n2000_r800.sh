#!/bin/bash

n=2000
r=800
#K=(825 850 875 900 925 950 975 1000 1050 1100 1150 1200 1250 1300 1400 1500 1600 1700 1800 1900 2000)
K=(805 810 815 820)

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

