function plotOptCar(t, x, u, y, alpha, params)
%%
n= length(t);

subplot(1, 2, 2)
plot(x(1, :), x(2, :))
grid on
title('x y')

subplot(4, 2 , 1)
plot(t(1:end-1), u(1, :))
grid on
title('steering angle')

subplot(4, 2, 3)
plot(t(1:end-1), u(2, :))
grid on
title('Acceleration')

subplot(4, 2, 5)
plot(t, x(3, :)/pi*180)
grid on
title('Car orientation')

subplot(4, 2, 7)
plot(t, x(4, :))
grid on
title('Car speed')
