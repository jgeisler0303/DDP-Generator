reset

if (!exists("datafile")) datafile='default.dat'
if (!exists("outfile")) outfile='default.eps'

set term postscript eps color 
set output outfile

set grid
set xlabel
set ylabel
set title "x y"
unset key
plot datafile using 1:2 title '' with line

set multiplot layout 2,2 rowsfirst
#set multiplot layout 4,1 title "Auto-layout of stacked plots\n" font ",12"

set title "Steering angle"
set grid
set xlabel "steps"
set ylabel 
unset key
#unset xtics
plot datafile using 5 title '' with line

set title "Acceleration"
set grid
set xlabel "steps"
set ylabel 
unset key
plot datafile using 6 title '' with line

set title "Car orientation"
set grid
set xlabel "steps"
set ylabel 
unset key
plot datafile using 3 title '' with line

set title "Car speed"
set grid
set xlabel "steps"
set ylabel 
unset key
plot datafile using 4 title '' with line

unset multiplot
