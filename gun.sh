#!/bin/bash

gcc  main.c functions.c functions_rsd.c nlbias_integrals.c pt_integrals.c pt_kernels.c tns_integrals.c tns_kernels.c cubature.c structures.c  -O3 -lm -fopenmp -o file_RPT.out

a2=0
b2=$a2

while [ $a2 -le 2 ]; do

a3=0
while [ $a3 -le 9 ] ; do

a4=1
b4=0

b3=`expr $a3 + 1`

if [ $b3 -eq 10 ]
then

b3=0
b2=`expr $b2 + 1`

fi

k=`expr $a2$a3$a4`

k2=`expr $b2$b3$b4`

sed -e "s/REAL2/$k2/g"  params_reference.ini > params_template_temp1.ini

sed -e "s/REAL/$k/g"  params_template_temp1.ini > params_template_temp$k.ini

sed -e "s/REAL/$k/g" script_reference.sh  > script_temp$k.sh

echo "$k -> $k2"


sbatch -H script_temp$k.sh

rm  params_template_temp1.ini


a3=`expr $a3 + 1`

done

if [ $a3 -eq 9 ]
then
a3=0
fi

a2=`expr $a2 + 1`

done

if [ $a2 -eq 9 ]
then
a2=0
fi
