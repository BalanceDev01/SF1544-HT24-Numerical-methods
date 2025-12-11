clc, clear, close all
%coefficients and functions
m1 = 475;
m2 = 53;
k1 = 5400;
k2 = 135000;
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

tspan = [0 (L/u)+1];
y0 = zeros(4,1); % y0 is v0
n1 = round((tspan(2)-tspan(1))/(5*1e-3)); %ger (ungefär) tiddsteg ∆t = 5*10^-3
n2 = round((tspan(2)-tspan(1))/(5*1e-4)); %ger (ungefär) tiddsteg ∆t = 5*10^-4 %För exakt tiddsteg, gör en ny Eulerquartercar och ändra steglängd h till 5*10^-4

[tv1, yv1] = Eulerquartercar(f, tspan, y0, n1); % Euler method med ∆t = 5e-3
[tv2, yv2] = Eulerquartercar(f, tspan, y0, n2); % Euler method med ∆t = 5e-4

z2_n1 = yv1(2,:);
z2_n2 = yv2(2,:);

options = odeset('RelTol',1e-6,'AbsTol',1e-7, 'Refine',1);
[t,v] = ode45(f,tspan,y0,options);

plot (t,v(:,2),'Color','b')
hold on;
plot(tv1, z2_n1,'Color', '#0072BD', 'LineWidth', 2);
hold on;
plot(tv2, z2_n2,'Color', '#4DBEEE', 'LineWidth', 2);
xlabel('Tid[s]')
ylabel('Höjd[cm]')
legend('z1(ode45)','z2(∆t = 5*10^-3)','z1(∆t = 5*10^-4')

% Vi ser här att det mindre tidssteget 10^-4 i gav betydligt närmare resultat till
% lösningen ode45 än 10^-3 i tidssteg och ode45 ger bästa lösningen då denna
% är en blandning av 4:e och 5:e ordningens noggrannhet.