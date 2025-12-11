clc, clear, close all
%coefficients and functions
m1 = 475;
m2 = 53;
k1 = 5400;
k2 = 135000;
k2 = 100 * k2; % för uppgift 3d), kommentera bort för uppgift 3c)
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
f = @(t,v) A*v + g(t); % f is dv/dt

alphas = [0.9 1 1.1 1.5]

% 3c) Testar vi att köra här med alpha(1) ser vi att systemet är stabilt.
% Ändrar vi till alpha(2) så stabiliserar inte systemet utan fortsätter
% oscillera. Tar vi alpha(3) eller alpha(4) så divergerar systemet då t ->
% inf, som väntat!

% delta_t = alphas(1) * 0.010781  % Denna för uppgift 3c
delta_t = 0.1 * 0.00011181 % Denna för uppgift 3d
tspan = [0 (L/u)+1];
y0 = zeros(4,1); % y0 is v0
n = round(tspan(2)/delta_t);
[tv,yv] = Eulerquartercar(f,tspan,y0,n); %Euler method
yv;

z1 = yv(1,:); %chassi vertical pos
z2 = yv(2,:); %wheel vertical pos
plot(tv, z1,'Color', 'r', 'LineWidth', 2);
hold on;
plot(tv, h(tv),'Color', 'black', 'LineWidth', 2);
hold on;
plot(tv, z2,'Color', 'b', 'LineWidth', 2);
hold off;
legend('z1', 'h(t)', 'z2');






