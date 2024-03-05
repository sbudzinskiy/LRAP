#!/bin/bash

n=2000
r=1600
#K=(1625 1650 1675 1700 1725 1750 1775 1800 1825 1850 1875 1900 1925 1950 1975 2000)
K=(1605 1610 1615 1620)

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

