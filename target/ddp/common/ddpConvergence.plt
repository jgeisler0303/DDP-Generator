reset

if (!exists("datafile")) datafile='ddp_iter.dat'
if (!exists("outfile")) outfile='ddp_iter.eps'

set term postscript eps color 
set output outfile

set grid
set xlabel "Iterations"
set ylabel "Cost"
set y2label "log10 of rest"
set title "Convergence"
set ytics nomirror
set xtics nomirror
set y2tics
set tics out
set logscale y2 10
set key bottom left
#set autoscale  y
#set autoscale y2
plot datafile using 1 title "cost" with line axes x1y1, \
    datafile using 2 title "lambda" with line axes x1y2, \
    datafile using 3 title "alpha" with line axes x1y2, \
    datafile using 4 title "gradient" with line axes x1y2, \
    datafile using 5 title "improvement" with line axes x1y2
    

reset
set title "Reduction ratio"
set grid
set xlabel "Iterations"
set ylabel 
plot datafile using 6 title "Reduction ratio" with line


