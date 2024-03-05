#!/bin/bash

n=1500
R=(1400 1200 1000 800 600 400 200 100 50)

p=5
q=1

executable="../../build/intel-release-mkl/bin/matrix_identity"

for r in "${R[@]}"
do
    filename="n${n}_r${r}_rsvd_p${p}_q${q}"
    params="${n} ${r} 2 ${p} ${q}"
    for seed in {1..5}
    do
        mkdir -p "./${filename}"
	echo "n ${n} r ${r} seed ${seed}"
	$executable $seed $params > "./${filename}/${seed}.log"
    done
done

