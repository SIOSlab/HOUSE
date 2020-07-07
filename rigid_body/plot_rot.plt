set terminal epslatex size 8in,10in

set xlabel 'Time (s)'

set logscale y
set yrange [1E-4:1E1]

set format y '$10^{%T}$'

set key outside
set key right

set style line 4 dashtype 3 lw 4 lc rgb '#0072bd' # blue
set style line 2 dashtype 1 lw 4 lc rgb '#d95319' # orange
set style line 3 dashtype 2 lw 4 lc rgb '#77ac30' # green
set style line 1 dashtype 4 lw 4 lc rgb '#7e2f8e' # purple

set output 'rot_error.tex'

set multiplot layout 3,1

set ylabel '$\omega_1$ (rad/s)'
plot 'dxus' u 1:(abs($2-$5)) title 'UKF Error' w lines ls 1, \
     'dxhm' u 1:(abs($2-$5)) title 'HOUSE Error' w lines ls 2, \
     'dxus' u 1:8 title 'UKF Standard Deviation' w lines ls 3, \
     'dxhm' u 1:8 title 'HOUSE Standard Deviation' w lines ls 4

set ylabel '$\omega_2$ (rad/s)'
plot 'dxus' u 1:(abs($3-$6)) title 'UKF Error' w lines ls 1, \
     'dxhm' u 1:(abs($3-$6)) title 'HOUSE Error' w lines ls 2, \
     'dxus' u 1:9 title 'UKF Standard Deviation' w lines ls 3, \
     'dxhm' u 1:9 title 'HOUSE Standard Deviation' w lines ls 4

set ylabel '$\omega_3$ (rad/s)'
plot 'dxus' u 1:(abs($4-$7)) title 'UKF Error' w lines ls 1, \
     'dxhm' u 1:(abs($4-$7)) title 'HOUSE Error' w lines ls 2, \
     'dxus' u 1:10 title 'UKF Standard Deviation' w lines ls 3, \
     'dxhm' u 1:10 title 'HOUSE Standard Deviation' w lines ls 4

unset multiplot
