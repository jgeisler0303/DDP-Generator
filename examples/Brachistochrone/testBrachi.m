setenv('MAXIMA', 'path to maxima');
addpath('../..');
make_iLQG('optDefBrachi', '-DDEBUG_BACKPASS=1 -DDEBUG_FORWARDPASS=1 -DFULL_DDP=0', 1)


%% parameters
p.g= 9.81;
p.yf= -4;
p.w= 100;
x0= [-eps];
Op.max_iter= 20;

clf
hold on
for n= [2 3 5 500]
    p.dx= 2*pi/n;

    x= linspace(0, 2*pi, n+1);
    u0= -ones(1,n);      % initial controls

    [success, y, u, cost]= iLQGBrachi(x0, u0, p, Op);
    
    plot(x, y, 'b')
end

a= 2;
phi_end= pi;
phi= linspace(0, phi_end, 1000);
x_true= a*(phi-sin(phi));
y_true= a*(cos(phi)-1);

plot(x_true, y_true, 'r')
grid on
hold off