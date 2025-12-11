%% ode45 lösning
close all;
tic
t_start = 0;
t_slut = 500;
a= 0.5;
x_0 = [1-a; 0; 0; sqrt((1+a)/(1-a))];

%x = lösningsvektorn

r = @(t,x) (x(1)^2 + x(2)^2) ^ (3/2);
F = @(t,x) [x(3); x(4); -x(1)/r(t,x); -x(2)/r(t,x)];

options = odeset('RelTol', 1e-15, 'AbsTol', 1e-18);
[t, sol] = ode45(F, [t_start, t_slut], x_0, options);

q1 = sol(:,1);
q2 = sol(:,2);
p1 = sol(:,3);
p2 = sol(:,4);


% Plotta banan
figure;
plot(q1, q2, 'b', 'LineWidth', 1.5);
xlabel('q1');
ylabel('q2');
title('Bana i (q1, q2)-planet (ode45)');
grid on;
axis equal;


N = length(t);
E = zeros(1, N);
for n = 1:N
    q = [q1(n) ; q2(n)];
    p = [p1(n) ; p2(n)];
    E(n) = (1/2) * norm(q)^2 - 1 / norm(p);
end


figure;
plot(t, E, 'r', 'LineWidth', 1.5);
xlabel('Tid');
ylabel('Total energi');
title('Energins förändring över tid');
grid on;
last_sol_ode = [p1(N), p2(N), q1(N), q2(N)];
fprintf('Lösningsvektorn vid t = %.0f för ode45 är: \n', t(N))
disp(last_sol_ode)
%fprintf('Lösningsvektorn vid t = %.0f är %.2f sekunder.\n', )

elapsedTime = toc;
fprintf('ode45 tog %.2f sekunder.\n', elapsedTime)

%% Mittpunnktsmetoden
% Parametrar
t_start = 0;
t_slut = 500; % 100
h = [0.5, 0.05, 0.005, 0.0005, 0.00005];  % Steglängd
felen = zeros(1,length(h));

for i = 1:length(h)
    t = t_start:h(i):t_slut;
    N = length(t);
    
    % Initialvärden
    q1 = zeros(1, N);
    q2 = zeros(1, N);
    p1 = zeros(1, N);
    p2 = zeros(1, N);
    
    a = 0.5;
    q1(1) = 1 - a;
    q2(1) = 0;
    p1(1) = 0;
    p2(1) = sqrt((1 + a) / (1 - a));
    
    % Implicit Mittpunktsmetod med Newtons metod
    %Nogrannhetsordning = 2
    for n = 1:N-1
        % Startgissning för Newtons metod (medelvärde av föregående värde)
        q_next = [q1(n); q2(n); p1(n); p2(n)];
        
        % Newtons metod
        tol = 1e-8;
        maxiter = 100; %200
        iter = 0;
        
        while iter < maxiter
            iter = iter + 1;
            
            % Medelpunkt
            q_half = (q_next + [q1(n); q2(n); p1(n); p2(n)]) / 2;
            r = q_half(1)^2 + q_half(2)^2;
            
            % Funktion att lösa: F(y) = y_n+1 - y_n - h*f(y_med)
            F = [q_next(1) - q1(n) - h(i) * q_half(3);
                 q_next(2) - q2(n) - h(i) * q_half(4);
                 q_next(3) - p1(n) - h(i) * (-q_half(1) / r^(3/2));
                 q_next(4) - p2(n) - h(i) * (-q_half(2) / r^(3/2))];
    
            % Jakobimatris (deriverad med avseende på q_next)
            J = [1, 0, -h(i)/2, 0;
                 0, 1, 0, -h(i)/2;
                 h(i)/2 * (-1/r^(3/2) + 3*q_half(1)^2/r^(5/2)), 3*h(i)/2 * q_half(1)*q_half(2)/r^(5/2), 1, 0;
                 3*h(i)/2 * q_half(1)*q_half(2)/r^(5/2), h(i)/2 * (-1/r^(3/2) + 3*q_half(2)^2/r^(5/2)), 0, 1];
    
            % Newton-steg: Lös J * s = -F
            s = -J \ F;
            q_next = q_next + s;
    
            % Kontrollera konvergens
            if norm(s, inf) < tol
                break;
            end
        end
        
        % Om Newtons metod inte konvergerar, ge en varning
        if iter == maxiter
            warning('Newton-metoden nådde max antal iterationer vid steg %d', n);
        end
    
        % Uppdatera lösningen
        q1(n+1) = q_next(1);
        q2(n+1) = q_next(2);
        p1(n+1) = q_next(3);
        p2(n+1) = q_next(4);
    end

    las_sol_mittpunkt = [p1(N), p2(N), q1(N), q2(N)];
    disp(las_sol_mittpunkt)
    fel = norm(last_sol_ode - las_sol_mittpunkt);
    felen(i) = fel;
end

disp('För mittpunktsmetoden:')
for i= 1:length(felen) - 1
        p = log(felen(i)/felen(i+1)) / log(h(i)/h(i+1));
        fprintf("Iteration %d: p = %.8f\n", i, p);
end


%% Symplektisk Euler
tic
% Parametrar
t_start = 0;
t_slut = 500; %100
h = [0.5, 0.05, 0.005, 0.0005, 0.00005, 0.00003, 0.00001];  % Steglängd
felen = zeros(1,length(h));
disp('sista lösningsvektorn (p1,p2,q1,q2) för steglängden är: ')
for i = 1: length(h)
    t = t_start:h(i):t_slut;
    N = length(t);
    
    % Initialvärden
    q1 = zeros(1, N);
    q2 = zeros(1, N);
    p1 = zeros(1, N);
    p2 = zeros(1, N);
    
    a = 0.5;
    q1(1) = 1 - a;
    q2(1) = 0;
    p1(1) = 0;
    p2(1) = sqrt((1 + a) / (1 - a));
    
    % Symplektisk Euler-iteration
    %Nogrannhetsordning = 1 (Pga framåt euler används)
    for n = 1:N-1
        % Uppdatera impulsen först (p)
        r = (q1(n)^2 + q2(n)^2)^(3/2);
        p1(n+1) = p1(n) - h(i) * (q1(n) / r);
        p2(n+1) = p2(n) - h(i) * (q2(n) / r);
    
        % Uppdatera positionen med nya impulsen
        q1(n+1) = q1(n) + h(i) * p1(n+1);
        q2(n+1) = q2(n) + h(i) * p2(n+1);
        % if tic > 1e10
        %     error('KRASH');
        % end 
    end
    las_sol_symplEuler = [p1(N), p2(N), q1(N), q2(N)];
    disp(las_sol_symplEuler)
    fel = norm(last_sol_ode - las_sol_symplEuler);
    felen(i) = fel;
end 

disp('För symplektisk Euler:')
p_vec = zeros(1,length(h));
for i= 1:length(felen) - 1
        p = log(felen(i)/felen(i+1)) / log(h(i)/h(i+1));
        fprintf("Iteration %d: p = %.8f\n", i, p);
        p_vec(i) = p;
end

avreage = sum(p_vec(2:5))/5;
fprintf('Medelvärdet från andra iterationen till sista är %.8f\n', avreage);

elapsedTime = toc;
fprintf('Symplektisk euler tog %.2f sekunder.\n\n', elapsedTime)