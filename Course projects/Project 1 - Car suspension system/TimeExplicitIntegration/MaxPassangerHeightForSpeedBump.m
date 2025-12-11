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
options = odeset('RelTol',1e-6,'AbsTol',1e-7, 'Refine',4);
[t,v] = ode45(f,tspan,y0,options);

[max_z1, i_z1] = max(v(:,1)); %max höjd av z1
[max_z2, i_z2] = max(v(:,2)); %max höjd av z2
disp(['Max height z1: ', num2str(max_z1)]);
disp(['Max height z2: ', num2str(max_z2)]);
disp(['Passagerarna upplever en förskjutning i y-led på ',num2str(100*max_z1), 'cm'])

plot (t,v(:,1),'Color','r')
hold on;
plot (t,v(:,2),'Color','b')
hold on;
xlabel('Tid [s]')
ylabel('Höjd [cm]')
legend('z1','z2');
