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

tspan = [0 0.05];
y0 = zeros(4,1);

% Referens lösning med ode45
rel_tol = 1e-9;
abs_tol = 1e-9;
options = odeset('RelTol', rel_tol, 'AbsTol', abs_tol);

[t_ref, y_ref] = ode45(@(t, y) A * y + g(t), tspan, y0, options);
z2_ref = y_ref(:, 2);

% konvergens studien
dt0 = 0.001; % initial tidssteg
alphas = [1, 1/2, 1/4, 1/8];
errors = zeros(size(alphas));

for i = 1:length(alphas)
    delta_t = alphas(i) * dt0; % steg storlek
    n = round(tspan(2) / delta_t); % antal steg
    [tv, yv] = implicit_trapezoid(A, g, tspan, y0, n); % lös med implicit trapets
    z2_trap = interp1(tv, yv(2, :), t_ref, 'linear'); % interpolera så vi får punkter vid samma tider
    errors(i) = max(abs(z2_trap - z2_ref));
end

fprintf('Steg storlekar:\n');
disp(alphas * dt0);
fprintf('Fel normer:\n');
disp(errors);

% Beräkna noggrannhetsordning
orders = zeros(1, length(errors) - 1);
for i = 1:length(errors) - 1
    orders(i) = log(errors(i) / errors(i+1)) / log((alphas(i) * dt0) / (alphas(i+1) * dt0));
end

fprintf('Beräknade noggrannhetsordningar:\n');
disp(orders);


% Vi ser att den numeriskt beräknade norgrannhetsordningen är 2. 
% Det stämmer bra med teorin.
