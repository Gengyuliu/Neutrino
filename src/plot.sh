#!/bin/bash
Nt=$1
Nz=$2
Nx=$3
N=$(($Nz*$Nx))
for t in `seq 0 100 $Nt`
do
    file=`printf "./rho_%.5d.txt" "$t"`
    gnuplot --persist << EOF
        unset label
        set title "At n = $t, z = -7.95" 
        set xlabel "x"
        set ylabel "rho_{ee}"
        set terminal png
        set output 'rho_$t.png'
        set contour
        set dgrid3d $Nz, $Nx
        set pm3d
        #set logscale z 10
        #set zrange[1e-16:1e-9]
        #splot '$file' u 3:4:(sqrt(column(5)**2+column(6)**2)) every ::130001::140000 w l notitle
        #splot '$file' u 3:4:5 every ::1::60000
        plot '$file' u 4:5 every ::4001::4200 t "RK4+FD4", '' u 4:6 every ::4001::4200 t "exact"
        #set terminal pop
        #set output
        replot
EOF
done

