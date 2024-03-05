#!/bin/bash

m=1000
n=1000
r=5
os=7.0
niter_rgd=1000

# 1 -- AP, 2 -- RP, 3 -- IP, 4 -- SP
method_ap=1
niter_ap=1
shift_ratio=0.9

# 1 -- CROSS, 2 -- CROSS-P
method_lr=1
niter_maxvol=10
thresh=1.05
k=30

executable="../../build/intel-release-mklseq/bin/matrix_cmplt"
params="${m} ${n} ${r} ${os} ${niter_rgd} ${method_ap} ${niter_ap}"
filename="os${os}_nap${niter_ap}"
if [ "$method_ap" = 1 ]; then
    filename="${filename}_ap"
elif [ "$method_ap" = 2 ]; then
    filename="${filename}_rp"
elif [ "$method_ap" = 3 ]; then
    filename="${filename}_ip_${shift_ratio}"
    params="${params} ${shift_ratio}"
else
    filename="${filename}_sp_${shift_ratio}"
    params="${params} ${shift_ratio}"
fi

params="${params} ${method_lr}"
if [ $method_lr = 1 ]; then
    filename="${filename}_cross"
    params="${params} ${niter_maxvol} ${thresh}"
else
    filename="${filename}_crossp_k${k}"
    params="${params} ${niter_maxvol} ${thresh} ${k} ${k}"
fi

for seed in {1..20}
do
    mkdir -p "./${filename}"
    echo "seed ${seed}"
    $executable $seed $params > "./${filename}/${seed}.log"
done
