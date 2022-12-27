function [ynew, ymod, M] = solver(t, y, dt)
% Runge-Kutta method 3rd order

[s1, ymod, M] = equ(t, y, dt);
s2 = equ(t + dt, ymod + s1*dt, dt);
s3 = equ(t + (1/2)*dt, ymod + (s1 + s2)/4*dt, dt);
s = (s1 + 4*s3 + s2)/6;

ynew = ymod + s*dt;
