reset

stats 'mat/two_grid_residual_evolution.dat' using 2 nooutput
y_min = STATS_min
y_max = STATS_max

set title "Evolution of the residual of the two-grid method"
set key inside center bottom vertical Right noreverse enhanced autotitle box lt black linewidth 1.000 dashtype solid

set logscale y
set grid xtics ytics
set style data lines
set xtics border out scale 1,0.5 nomirror
set ytics log
set ytics border out scale 1,0.5 nomirror

set yrange [y_min/10:y_max*10]
set xlabel "Iteration"
set ylabel "Residual"
plot "mat/two_grid_residual_evolution.dat" with lines
replot