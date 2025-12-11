clear; clc; close all;

N=200;       % antal intervall
T=1;         % sluttid
dx=1/N;    % steglängd i rummet
dt=dx/2.0; % tidssteg, tänk på stabilitetsvillkoren
M=round(T/dt); % antal tidsteg
c=2;

                  % allokering av minne
u=zeros(N-1,M+1); % u(n,m): lösningens värde vid tid (m-1)*dt i position n*dx
p=zeros(N-1,M+1); % p=u'

A=zeros(N-1,N-1);  % Au differensapproximation av d^2 u/dx^2
x = dx*(1:N-1)';   % x(n) är n*dx
E = zeros(1,M+1);  % För att beräkna energin i varje tidssteg.

%Skapa matrisen A
mittenterm = -2/dx^2 * ones(N-1,1);
A = diag(mittenterm);
sidotermer = 1/dx^2 * ones(N-2,1);
A = A + diag(sidotermer,1) + diag(sidotermer,-1);
%disp(A)

%Sätt begynnelsedata för u och p.
g = @(x) exp(-200*(x-0.5).^2);
u(:,1) = g(x);
p(:,1) = zeros(size(x));

%Räkna ut energin E vid tiden 0.
%E=zeros(1,M);
E(1) = 0;
Hamiltonian = @(p,q) (1/2)*(norm(p)^2 - dot(q,c^2*A*q)); %c^2*A*q


nframe=M+1;   % kommando för film
mov(1:nframe)= struct('cdata', [], 'colormap', []);
figure;
subplot(1,2,1)
plot(x, u(:, 1), 'b', 'LineWidth', 1); %Plot vid tiden t=0.
axis([0 1 -1 1]);
set(gca, 'nextplot', 'replacechildren');
title('Lösning vid tiden t');
drawnow
mov(1)=getframe(gcf); %Första frame i filmen.


for m = 1:M        % tidsstegning med symplektiska Euler
    u(:, m+1) = u(:,m) + dt * p(:,m);
    p(:, m+1) = p(:,m) + dt * c^2 * A * u(:, m+1); 
    u_nu = u(:, m+1);
    p_nu = p(:, m+1);

    X = [0; x; 1]; U = [0; u(:, m); 0];
    
    subplot(1,2,1)
    plot(X, U, 'b', 'LineWidth', 1);
    hold on;
    %Plotta även lösningen från d'Alemberts formel

    t = dt*(m-1);
    d_Lambert = (1/2)*(g(x - c*t) + g(x + c*t));
    subplot(1,2,1)
    plot(X, [0;d_Lambert;0], 'r', 'LineWidth', 1);
    hold on;

    text(0.05,-0.8, sprintf('t=%.2f', m * dt));
    set(gca, 'nextplot', 'replacechildren');
    drawnow
    mov(m+1) = getframe(gcf);
    
    %Räkna ut energin av den numeriska lösningen vid detta tidsstag
    E(m+1) = Hamiltonian(p_nu, u_nu);
    disp(E(m+1)) %visar energins värden, förhoppningsvis typ konstant
    subplot(1,2,2);
    t_vec = dt * (0:m);
    plot(t_vec,E(1,1:m+1), 'r', 'LineWidth', 1);
    hold on;
end

movie(mov, 1, 30);
subplot(1,2,2)
t_axis=dt*(0:M);
plot(t_axis, E(1), 'b', 'LineWidth', 1);
axis([0 2 6.5e3 7.5e3]); 
set(gca, 'nextplot', 'replacechildren');
title('Energin');
drawnow
mov(1) = getframe(gcf); %Första frame i filmen.

