clc, clear, close all
%coefficients and functions
m1 = 475;
m2 = 53;
k1 = 5400;
k2 = 135000 * 100;
c1 = 310;
c2 = 1200;
u = 65/3.6;
H = 0.24;
L = 1;
A = [0 0 1 0;
    0 0 0 1;
    -(k1/m1) k1/m1 -c1/m1 c1/m1;
    k1/m2 -(k1+k2)/m2 c1/m2 -(c1+c2)/m2];
h = @(t) (t <= L/u) .* (H/2 * (1 - cos(2 * pi * u * t / L))) + ...
         (t > L/u) .* 0;
dhdt = @(t) (t <= L/u) .* (H * pi * u / L * sin(2 * pi * u * t / L)) + ...
            (t > L/u) .* 0;
g = @(t) [0;
    0;
    0;
    (k2/m2) * h(t) + (c2/m2) * dhdt(t)];

alphas = [1 10 100];
delta_t = alphas(1) * 0.00011181;
tspan = [0 (L/u)+1];
y0 = zeros(4,1); % y0 is v0
n = round(tspan(2)/delta_t);
[tv,yv] = implicit_trapezoid(A,g,tspan,y0,n);

z1 = yv(1,:); %chassi vertical pos
z2 = yv(2,:); %wheel vertical pos
plot(tv, z1,'Color', 'r', 'LineWidth', 2);
hold on;
plot(tv, h(tv),'Color', 'black', 'LineWidth', 2);
hold on;
plot(tv, z2,'Color', 'b', 'LineWidth', 2);
hold off;
legend('z1', 'h(t)', 'z2');

% Lösningarna förblir stabila även med mycket längre tidssteg vilket är vad
% teorin säger.
