#!/bin/bash

eps="1e-2"

# 1 -- AP, 2 -- RP, 3 -- IP, 4 -- SP
method_ap=1
niter_ap=10
shift_ratio=0.5

#1 -- TTSVD, 2 -- TTCROSS
method_lr=1
p=0
q=0
niter_maxvol=20
thresh=1.1
k=30

executable="../../build/intel-release-mklseq/bin/tensor_gaussian"
params="${eps} ${method_ap} ${niter_ap}"
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
    filename="${filename}_ttsvd"
elif [ $method_lr = 2 ]; then
    filename="${filename}_ttcross"
    params="${params} ${niter_maxvol} ${thresh}"
fi

for seed in {1..1}
do
    mkdir -p "./${filename}"
    echo "seed ${seed}"
    $executable $seed $params > "./${filename}/${seed}.log"
done
