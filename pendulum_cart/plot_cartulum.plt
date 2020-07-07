set terminal epslatex size 6in,7.5in

set xlabel 'Time (s)'

set logscale y
set format y '$10^{%T}$

set key outside
set key right top

set style line 1 dashtype 1 lc 'blue'
set style line 2 dashtype 1 lc 'red'
set style line 3 dashtype 2 lc 'blue'
set style line 4 dashtype 2 lc 'red'

#-------------------------------------------------------------------------------

set output 'out/cartulum_error.tex'

set multiplot layout 4,1

set ylabel '$x_c$ (m)'
plot 'out/dxus' u 1:(abs($2-$6)) title 'UKF $\epsilon$' w lines ls 1, \
     'out/dxhm' u 1:(abs($2-$6)) title 'HOUSE $\epsilon$' w lines ls 2, \
     'out/dxus' u 1:10 title 'UKF $\sigma$' w lines ls 3, \
     'out/dxhm' u 1:10 title 'HOUSE $\sigma$' w lines ls 4

set ylabel '$\theta$ (rad)'
plot 'out/dxus' u 1:(abs(2*asin(sin(($3-$7)/2)))) title 'UKF $\epsilon$' w lines ls 1, \
     'out/dxhm' u 1:(abs(2*asin(sin(($3-$7)/2)))) title 'HOUSE $\epsilon$' w lines ls 2, \
     'out/dxus' u 1:11 title 'UKF $\sigma$' w lines ls 3, \
     'out/dxhm' u 1:11 title 'HOUSE $\sigma$' w lines ls 4

set ylabel '$\dot{x}_c$ (m/s)'
plot 'out/dxus' u 1:(abs($4-$8)) title 'UKF $\epsilon$' w lines ls 1, \
     'out/dxhm' u 1:(abs($4-$8)) title 'HOUSE $\epsilon$' w lines ls 2, \
     'out/dxus' u 1:12 title 'UKF $\sigma$' w lines ls 3, \
     'out/dxhm' u 1:12 title 'HOUSE $\sigma$' w lines ls 4

set ylabel '$\dot{\theta}$ (rad/s)'
plot 'out/dxus' u 1:(abs($5-$9)) title 'UKF $\epsilon$' w lines ls 1, \
     'out/dxhm' u 1:(abs($5-$9)) title 'HOUSE $\epsilon$' w lines ls 2, \
     'out/dxus' u 1:13 title 'UKF $\sigma$' w lines ls 3, \
     'out/dxhm' u 1:13 title 'HOUSE $\sigma$' w lines ls 4

unset multiplot

#-------------------------------------------------------------------------------
