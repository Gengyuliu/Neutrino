#!/bin/bash
theta=33.3
mu=1.0
f_0=1.0
alpha=0.92
Nt=400
CFL=0.2
dim=1
Nz=200
#dz=12
z1=-10.0
z2=10.0
Nvz=10
Nx=200
#dx=24
x1=-10.0
x2=10.0
Nvx=10
init=2      #another Gaussian
A=`echo "1/10000000" | bc -l`
#z_0=$(echo "$dz*$M_z/2.0" | bc -l)
z_0=0.0
#sigma=$(echo "$dz*$M_z/20.0" | bc -l)
sigma=2.0
#echo $z_0
time ./main --theta $theta --mu $mu --f_0 $f_0 --alpha $alpha --Nt $Nt --CFL $CFL --dim $dim --Nz $Nz --z1 $z1 --z2 $z2 --Nvz $Nvz  --A $A --z_0 $z_0 --sigma $sigma
#time ./main --theta $theta --mu $mu --f_0 $f_0 --alpha $alpha --Nt $Nt --CFL $CFL --dim $dim --Nz $Nz --z1 $z1 --z2 $z2 --Nvz $Nvz --Nx $Nx --x1 $x1 --x2 $x2 --Nvx $Nvx --A $A --z_0 $z_0 --sigma $sigma
# gdb --args ./main --theta 33.3 --Nt 10 --CFL 0.2 --dim 2 --Nz 200 --z1 -10 --z2 10 --Nvz 2 --Nx 200 --x1 -10 --x2 10 --Nvx 2 --z_0 0.0 --sigma 2.0



#./plot.sh $Nt $Nz $Nx
