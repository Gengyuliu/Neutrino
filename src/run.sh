#!/bin/bash
theta=0
mu=1.0
f_0=1.0
alpha=0.92
Nt=100
dt=0.5
dim=2
Nz=100
#dz=12
z1=-600.0
z2=600.0
Nvz=10
Nx=100
#dx=24
x1=-600.0
x2=600.0
Nvx=10
init=2      #another Gaussian
A=`echo "1/10000000" | bc -l`
#z_0=$(echo "$dz*$M_z/2.0" | bc -l)
z_0=0.0
#sigma=$(echo "$dz*$M_z/20.0" | bc -l)
sigma=5.0
#echo $z_0
#time ./main --theta $theta --mu $mu --f_0 $f_0 --alpha $alpha --Nt $Nt --dt $dt --dim $dim --Nz $Nz --z1 $z1 --z2 $z2 --Nvz $Nvz  --A $A --z_0 $z_0 --sigma $sigma
time ./main --theta $theta --mu $mu --f_0 $f_0 --alpha $alpha --Nt $Nt --dt $dt --dim $dim --Nz $Nz --z1 $z1 --z2 $z2 --Nvz $Nvz --Nx $Nx --x1 $x1 --x2 $x2 --Nvx $Nvx --A $A --z_0 $z_0 --sigma $sigma
./plot.sh $Nt $Nz $Nx
