
h = 0.001;  % Tidssteg
T = 100;    % Total tid
N = T/h;  % Antal steg
a = 0.5; % excentriciteten
beta= 0; % 1 för symplektisk euler då ska även alpha = 1/2

alpha = 0; %0 för euler bakåt, 1 för euler framåt, 1/2 för mittpunksmetoden


%Initialvärden (position och hastighet)
q1_0 = 1-a;
q2_0 = 0.0;
p1_0 = 0.0;
p2_0 = (1+a/1-a)^0.5;

%lagring av svar
q1_vals = zeros(1, N);
q2_vals = zeros(1, N);
p1_vals = zeros(1, N);
p2_vals = zeros(1, N);

% Initialisering
q1_vals(1) = q1_0;
q2_vals(1) = q2_0;
p1_vals(1) = p1_0;
p2_vals(1) = p2_0;

for n = 1:N-1

    r2 = q1_vals(n)^2 + q2_vals(n)^2;
    

    acc1 = -q1_vals(n+beta) / r2^(3/2);
    acc2 = -q2_vals(n+beta) / r2^(3/2);
    
    p1_vals(n+1) = p1_vals(n) + h * acc1;  
    p2_vals(n+1) = p2_vals(n) + h * acc2;  

    q1_vals(n+1) = (1-alpha) * (q1_vals(n) + h * p1_vals(n+1)) +  (alpha) * (q1_vals(n) + h * p1_vals(n));
    q2_vals(n+1) = (1-alpha) * (q2_vals(n) + h * p2_vals(n+1)) + (alpha) * (q2_vals(n) + h * p2_vals(n));

end

figure;
plot(q1_vals, q2_vals, 'red');
title('Bana för Keplerproblemet (mittpunktsmetod)');
xlabel('q1');
ylabel('q2');
axis equal;
