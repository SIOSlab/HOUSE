set terminal epslatex size 6in,7.5in

set xlabel 'Time (s)'

set logscale y
set format y '$10^{%T}$'

set key outside
set key right top

set style line 1 dashtype 1 lc 'blue'
set style line 2 dashtype 1 lc 'red'
set style line 3 dashtype 2 lc 'blue'
set style line 4 dashtype 2 lc 'red'

#-------------------------------------------------------------------------------

set output 'plots/proj_pos_error_hf.tex'

set multiplot title layout 3,1

set ylabel '$x$ (m)'
plot 'out/dxus' u 1:(abs($2-$8)) title 'UKF $\epsilon$' w lines ls 1, \
     'out/dxhm' u 1:(abs($2-$8)) title 'HOUSE $\epsilon$' w lines ls 2, \
     'out/dxus' u 1:14 title 'UKF $\sigma$' w lines ls 3, \
     'out/dxhm' u 1:14 title 'HOUSE $\sigma$' w lines ls 4

set ylabel '$y$ (m)'
plot 'out/dxus' u 1:(abs($3-$9)) title 'UKF $\epsilon$' w lines ls 1, \
     'out/dxhm' u 1:(abs($3-$9)) title 'HOUSE $\epsilon$' w lines ls 2, \
     'out/dxus' u 1:15 title 'UKF $\sigma$' w lines ls 3, \
     'out/dxhm' u 1:15 title 'HOUSE $\sigma$' w lines ls 4

set ylabel '$z$ (m)'
plot 'out/dxus' u 1:(abs($4-$10)) title 'UKF $\epsilon$' w lines ls 1, \
     'out/dxhm' u 1:(abs($4-$10)) title 'HOUSE $\epsilon$' w lines ls 2, \
     'out/dxus' u 1:16 title 'UKF $\sigma$' w lines ls 3, \
     'out/dxhm' u 1:16 title 'HOUSE $\sigma$' w lines ls 4

unset multiplot

#-------------------------------------------------------------------------------

set output 'plots/proj_vel_error_hf.tex'

set multiplot layout 3,1

set ylabel '$\dot{x}$ (m/s)'
plot 'out/dxus' u 1:(abs($5-$11)) title 'UKF $\epsilon$' w lines ls 1, \
     'out/dxhm' u 1:(abs($5-$11)) title 'HOUSE $\epsilon$' w lines ls 2, \
     'out/dxus' u 1:17 title 'UKF $\sigma$' w lines ls 3, \
     'out/dxhm' u 1:17 title 'HOUSE $\sigma$' w lines ls 4

set ylabel '$\dot{y}$ (m/s)'
plot 'out/dxus' u 1:(abs($6-$12)) title 'UKF $\epsilon$' w lines ls 1, \
     'out/dxhm' u 1:(abs($6-$12)) title 'HOUSE $\epsilon$' w lines ls 2, \
     'out/dxus' u 1:18 title 'UKF $\sigma$' w lines ls 3, \
     'out/dxhm' u 1:18 title 'HOUSE $\sigma$' w lines ls 4

set ylabel '$\dot{z}$ (m/s)'
plot 'out/dxus' u 1:(abs($7-$13)) title 'UKF $\epsilon$' w lines ls 1, \
     'out/dxhm' u 1:(abs($7-$13)) title 'HOUSE $\epsilon$' w lines ls 2, \
     'out/dxus' u 1:19 title 'UKF $\sigma$' w lines ls 3, \
     'out/dxhm' u 1:19 title 'HOUSE $\sigma$' w lines ls 4

unset multiplot
unset logscale

#-------------------------------------------------------------------------------

unset logscale y
unset format y

set xrange [0:12500]
set yrange [-12500:5000]

set xtics 0,5000
set ytics -10000,5000
set ztics 500

set terminal epslatex size 6in,6in

set output 'plots/proj_pos_true.tex'

set ticslevel 0

set xlabel '$x$ (m)$
set ylabel '$y$ (m)$
set zlabel '$z$ (m)$

splot 'out/dxus' u 2:3:4 notitle w lines ls 1

#-------------------------------------------------------------------------------


