#!/bin/bash
theta=33.3
mu=1.0
f_0=1.0
alpha=0.92
Nt=20
CFL=0.1
dim=2
Nz=3000
#dz=12
z1=-600.0
z2=600.0
Nvz=2
Nx=1000
#dx=24
x1=-100.0
x2=100.0
Nvx=2
init=2      #another Gaussian
A=`echo "1/10000000" | bc -l`
#z_0=$(echo "$dz*$M_z/2.0" | bc -l)
z_0=0.0
#sigma=$(echo "$dz*$M_z/20.0" | bc -l)
sigma=100.0
#echo $z_0
#time ./main --theta $theta --mu $mu --f_0 $f_0 --alpha $alpha --Nt $Nt --CFL $CFL --dim $dim --Nz $Nz --z1 $z1 --z2 $z2 --Nvz $Nvz  --A $A --z_0 $z_0 --sigma $sigma
time ./main --theta $theta --mu $mu --f_0 $f_0 --alpha $alpha --Nt $Nt --CFL $CFL --dim $dim --Nz $Nz --z1 $z1 --z2 $z2 --Nvz $Nvz --Nx $Nx --x1 $x1 --x2 $x2 --Nvx $Nvx --A $A --z_0 $z_0 --sigma $sigma
#./plot.sh $Nt $Nz $Nx
