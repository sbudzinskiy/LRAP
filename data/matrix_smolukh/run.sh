#!/bin/bash

r=10

# 1 -- AP, 2 -- RP, 3 -- IP, 4 -- SP
method_ap=4
niter_ap=1000
shift_ratio=0.5

#1 -- SVD, 2 -- RSVD, 3 -- CROSS, 4 -- CROSS-P
method_lr=2
p=2
q=0
niter_maxvol=10
thresh=1.05
k=30

executable="../../build/intel-release-mklseq/bin/matrix_smolukh"
params="${r} ${method_ap} ${niter_ap}"
if [ "$method_ap" = 1 ]; then
    filename="ap"
elif [ "$method_ap" = 2 ]; then
    filename="rp"
elif [ "$method_ap" = 3 ]; then
    filename="ip_${shift_ratio}"
    params="${params} ${shift_ratio}"
else
    filename="sp_${shift_ratio}"
    params="${params} ${shift_ratio}"
fi

params="${params} ${method_lr}"
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

for seed in {1..10}
do
    mkdir -p "./${filename}"
    echo "seed ${seed}"
    $executable $seed $params > "./${filename}/${seed}.log"
done
