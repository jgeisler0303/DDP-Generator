% setenv('MAXIMA', 'path to maxima');
% addpath('../..');
% make_target('optDefBrachi_hli', 'ddp/plain', '-DDEBUG_BACKPASS=1 -DDEBUG_FORWARDPASS=1 -DFULL_DDP=0', 1)


%% parameters
p.g= 9.81;
p.eps= 0;

x0= [-eps];
Op= [];
Op.max_iter= 20;
Op.w_pen_init_l= 40;
Op.w_pen_init_f= 1e-5;
Op.w_pen_max_f= 1;
Op.w_pen_fact2= 1;

clf
hold on
n= 500;
p.dx= 2*pi/n;
p.ymin= [linspace(-1, -5, n), -4];

x= linspace(0, 2*pi, n+1);
u0= -ones(1,n);      % initial controls

[success, y, u, cost]= ddpBrachi_hli(x0, u0, p, Op);

plot(x, y, 'b', x, p.ymin, 'r')

a= 2;
phi_end= pi;
phi= linspace(0, phi_end, 1000);
x_true= a*(phi-sin(phi));
y_true= a*(cos(phi)-1);

plot(x_true, y_true, 'g')
grid on
hold off