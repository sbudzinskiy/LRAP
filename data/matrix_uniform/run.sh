#!/bin/bash

n=250
r=200

# 1 -- SVD, 2 -- RSVD, 3 -- CROSS, 4 -- CROSS-P
method_lr=2
p=5
q=1
niter_maxvol=10
thresh=1.05
k=30

executable="../../build/intel-release-mklseq/bin/matrix_hadamard"
filename="n${n}_r${r}"
params="${n} ${r} ${method_lr}"
if [ $method_lr = 1 ]; then
    filename="${filename}_svd"
elif [ $method_lr = 2 ]; then
    filename="${filename}_rsvd_p${p}_q${q}"
    params="${params} ${p} ${q}"
elif [ $method_lr = 3 ]; then
    filename="${filename}_cross"
    params="${params} ${niter_maxvol} ${thresh}"
else
    filename="${filename}_crossp_k${k}"
    params="${params} ${niter_maxvol} ${thresh} ${k} ${k}"
fi

for seed in {1..1}
do
    mkdir -p "./${filename}"
    echo "seed ${seed}"
    $executable $seed $params > "./${filename}/${seed}.log"
done
