reset

if (!exists("datafile")) datafile='default.dat'
if (!exists("outfile")) outfile='default.eps'

set parametric
set samples 100
a= 2
#phi_end= pi
#phi= linspace(0, phi_end, 1000)
x_true(t)= a*(t-sin(t))
y_true(t)= a*(cos(t)-1)
set trange [0:pi]

constr(t)= -1 + -4*2*t/2/pi
constr_x(t)= 2*t

set term postscript eps color 
set output outfile

set xtics nomirror
set x2tics
set x2tics out
set x2range [0:2*pi]

set grid
set xlabel "steps"
set x2label "true solution x"
set ylabel 
plot datafile using 1 title 'ddp' with line, x_true(t),y_true(t)  title 'true solution' with line axes x2y1, constr_x(t),constr(t) title 'constraint' with line axes x2y1

