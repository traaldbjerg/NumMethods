reset

stats 'mat/time_evolution.dat' using 2 nooutput
y_min = STATS_min
y_max = STATS_max

set title "Evolution of the solution time in function of the number of unknowns of the problem, PCG with V-cycle preconditioner"
set key inside center bottom vertical Right noreverse enhanced autotitle box lt black linewidth 1.000 dashtype solid

#set logscale y
set grid xtics ytics
set style data lines
set xtics border out scale 1,0.5 nomirror
set ytics log
set ytics border out scale 1,0.5 nomirror

set yrange [-2 :y_max + 2]
set xlabel "Unknowns"
set ylabel "Solution time (s)"
plot "mat/time_evolution.dat" with lines
replot
replot