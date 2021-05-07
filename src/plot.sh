#!/bin/bash
Nt=$1
Nz=$2
Nx=$3
N=$(($Nz*$Nx))
for t in `seq 0 10 $Nt`
do
    file=`printf "./P_%.5d.txt" "$t"`
    gnuplot --persist << EOF
        unset label
        set xlabel "x"
        set ylabel "z"
        set terminal png
        set output 'P_$t.png'
        set contour
        set dgrid3d $Nz, $Nx
        set pm3d
        #set logscale z 10
        #set zrange[1e-16:1e-9]
        splot '$file' u 3:4:(sqrt(column(5)**2+column(6)**2)) every ::130001::140000 w l notitle
        #set terminal pop
        #set output
        #replot
EOF
done

