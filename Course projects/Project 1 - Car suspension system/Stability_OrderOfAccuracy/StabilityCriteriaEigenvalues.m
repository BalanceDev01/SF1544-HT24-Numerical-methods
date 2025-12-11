% 3a)
% Stabilitets villkoret ges av |1 + del_t * lambda| < 1 (enligt botten 
% på s2 i absolutstabilitet.pdf), där detta skall gälla för alla egenvärden 
% lambda av "systemmatrisen" A för systemet. Alltså måste vi ta det minsta
% tidssteget (del_t = h) beräknat på det här viset för alla egenvärden av A.

% 3b)
% vi fortsätter genom att beräkna egenvärden för A med samma parametrar som
% för uppgift 2. Sedan löser vi olikheten för del_t för varje egenvärde av
% A, och tar det minsta. Då får vi det största tillåtna tiddsteget.

clc, clear, close all

%parametrar
m1 = 475;
m2 = 53;
k1 = 5400;
k2 = 135000;
% k2 = 100 * k2; %avkommentera för fråga 3d)
c1 = 310;
c2 = 1200;
u = 65/3.6;
H = 0.24;
L = 1;
A = [0 0 1 0;
    0 0 0 1;
    -(k1/m1) k1/m1 -c1/m1 c1/m1;
    k1/m2 -(k1+k2)/m2 c1/m2 -(c1+c2)/m2];

lambdas = eig(A);
dt = inf;

for k = 1:length(lambdas)
    lambda = lambdas(k)
    alpha = real(lambda)
    beta = imag(lambda)
    a = alpha^2 + beta^2
    b = 2 * alpha
    dt_k = -b/a
    dt = min(dt, dt_k)
end

disp('Egenvärden för A:');
disp(lambdas);
disp(['Max tillåtna tidssteg: ', num2str(dt)]);

% 3b forts.) Eftersom tidssteget blir 0.010781 så uppfyller 0.005 från uppgift 2c 
% kravet!

% 3d) Om vi multiplicerar ursprungliga k2 med 100 så får vi max tidssteg
% till 0.00011181.

