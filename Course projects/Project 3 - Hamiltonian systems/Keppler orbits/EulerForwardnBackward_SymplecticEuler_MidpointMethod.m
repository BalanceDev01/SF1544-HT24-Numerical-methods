a=0.5; %excentricitet
tspan=[0 100];
h=0.0005;
n=200000;%=(tspan(2)-tspan(1))/h; 
tv=(tspan(1)+h*(0:n)); %tv, row vector. 
y0=[1-a; 0; 0;sqrt((1+a)/(1-a))];
yv=zeros(length(y0),n+1); % matrix with row vectors
yv(:,1) = y0;

%Eulerfram
for ii = 1:n 
    yv(:,ii+1) = [yv(1,ii) + h * yv(3,ii);
    yv(2,ii) + h * yv(4,ii);
    yv(3,ii) + h * (-yv(1,ii) / (yv(1,ii)^2 + yv(2,ii)^2)^(3/2));
    yv(4,ii) + h * (-yv(2,ii) / (yv(1,ii)^2 + yv(2,ii)^2)^(3/2))];
end

E = zeros(1,n+1); %energi 
for ii = 1:n+1
    r = sqrt(yv(1,ii)^2 + yv(2,ii)^2);
    E(ii) = 0.5 * (yv(3,ii)^2 + yv(4,ii)^2) - 1 / r;
end

%plottar omloppsbana
figure;
plot(yv(1,:), yv(2,:), 'b');
xlabel('q1');
ylabel('q2');
title('Omloppsbana, Euler fram');
grid on;
axis equal;

%plottar energi
figure;
plot(tv, E, 'r');
xlabel('Tid');
ylabel('Total energi');
title('Energins förändring över tid, Euler fram');
grid on;

%% 
a=0.5; %excentricitet
tspan=[0 100];
h=0.0005;
n=200000;%=(tspan(2)-tspan(1))/h; 
tv=(tspan(1)+h*(0:n)); %tv, row vector. 
y0=[1-a; 0; 0;sqrt((1+a)/(1-a))];
yv=zeros(length(y0),n+1); % matrix with row vectors
yv(:,1) = y0;

%Eulerbak
for ii = 1:n

    y_i1 = [yv(1,ii); yv(2,ii); yv(3,ii); yv(4,ii)]; %startgissning(av nästa steg)
    
    diffy=1; iter=0; maxiter=100; tol=10^(-8);
    
    while iter < maxiter
        iter = iter + 1;
        
        r = sqrt(y_i1(1)^2 + y_i1(2)^2); %y_i1 är nästa steg vi vill lösa för
        F = [y_i1(1) - yv(1,ii) - h * y_i1(3);
             y_i1(2) - yv(2,ii) - h * y_i1(4);
             y_i1(3) - yv(3,ii) - h * (-y_i1(1) / r^3);
             y_i1(4) - yv(4,ii) - h * (-y_i1(2) / r^3)];

        DF = [1, 0, -h, 0;
            0, 1, 0, -h;
            h*((1/(r^3)) - 3*(y_i1(1)^2)/(r^5)), h*(-3)*y_i1(1)*y_i1(2)/(r^5), 1, 0;
            h*(-3)*y_i1(1)*y_i1(2)/(r^5), h*((1/(r^3)) - 3*(y_i1(2)^2)/(r^5)), 0, 1];
        s = -DF \ F;
        y_i1 = y_i1 + s;

        if iter == maxiter
            disp(['max iterationer vid steg', num2str(n)]);
        end

        if norm(s, inf) < tol %bryter om newton ej konv
            break;
        end
    end

    %uppdaterar
    yv(:,ii+1) = [y_i1(1);y_i1(2); y_i1(3); y_i1(4)];
end

%plottar bana
figure;
plot(yv(1,:), yv(2,:), 'b');
xlabel('q1');
ylabel('q2');
title('Omloppsbana, Euler bak');
grid on;
axis equal;

E = zeros(1,n+1); %energi 
for ii = 1:n+1
    r = sqrt(yv(1,ii)^2 + yv(2,ii)^2);
    E(ii) = 0.5 * (yv(3,ii)^2 + yv(4,ii)^2) - 1 / r;
end

figure;
plot(tv, E, 'r');
xlabel('Tid');
ylabel('Total energi');
title('Energins förändring över tid, Euler bak');
grid on;

%% 
a=0.5; %excentricitet
tspan=[0 100];
h=0.00005;
n=2000000;%=(tspan(2)-tspan(1))/h; 
tv=(tspan(1)+h*(0:n)); %tv, row vector. 
y0=[1-a; 0; 0;sqrt((1+a)/(1-a))];
yv=zeros(length(y0),n+1); % matrix with row vectors
yv(:,1) = y0;

%Mittpunkt
for ii = 1:n

    y_i1 = [yv(1,ii); yv(2,ii); yv(3,ii); yv(4,ii)]; %startgissning(av nästa steg)
    
    diffy=1; iter=0; maxiter=100; tol=10^(-8);
    
    while iter < maxiter
        iter = iter + 1;

        y_i = [yv(1,ii); yv(2,ii); yv(3,ii); yv(4,ii)];
        r = sqrt((y_i1(1)+y_i(1))^2 + (y_i1(2)+y_i(2))^2); %y_i1 är nästa steg vi vill lösa för
        F = [y_i1(1) - yv(1,ii) - h * 0.5*(y_i1(3)+y_i(3));
             y_i1(2) - yv(2,ii) - h * 0.5*(y_i1(4)+y_i(4));
             y_i1(3) - yv(3,ii) + h * 4*(y_i1(1)+y_i(1))/r^(3);
             y_i1(4) - yv(4,ii) + h * 4*(y_i1(2)+y_i(2))/r^(3)];

        DF = [1, 0, -h/2, 0;
            0, 1, 0, -h/2;
            h*4*((1/r^3) + 2*(-3)*y_i1(1)*(y_i1(1)+y_i(1))/r^5), h*4*-3*(y_i1(1)+y_i(1))*(y_i1(2)+y_i(2))/r^5, 1, 0;
            h*4*-3*(y_i1(1)+y_i(1))*(y_i1(2)+y_i(2))/r^5, h*4*((1/r^3) + 2*(-3)*y_i1(2)*(y_i1(2)+y_i(2))/r^5), 0, 1];
        s = -DF \ F;
        y_i1 = y_i1 + s; %uppdaterar newtoniteration för lösn av y_i1

        if iter == maxiter
            disp(['max iterationer vid steg', num2str(n)]);
        end

        if norm(s, inf) < tol %bryter om newton ej konv
            break;
        end
    end

    %uppdaterar steg y_i (
    yv(:,ii+1) = [y_i1(1);y_i1(2); y_i1(3); y_i1(4)];
end

%plottar bana
figure;
plot(yv(1,:), yv(2,:), 'b');
xlabel('q1');
ylabel('q2');
title('Omloppsbana, mittpunktsmetod');
grid on;
axis equal;

E = zeros(1,n+1); %energi 
for ii = 1:n+1
    r = sqrt(yv(1,ii)^2 + yv(2,ii)^2);
    E(ii) = 0.5 * (yv(3,ii)^2 + yv(4,ii)^2) - 1 / r;
end

figure;
plot(tv, E, 'r');
xlabel('Tid');
ylabel('Total energi');
title('Energins förändring över tid, mittpunktsmetod');
grid on;
% ylim([-0.51 -0.49])
%%
a=0.5; %excentricitet
tspan=[0 100];
h=0.0005;
n=200000;%=(tspan(2)-tspan(1))/h; 
tv=(tspan(1)+h*(0:n)); %tv, row vector. 
y0=[1-a; 0; 0;sqrt((1+a)/(1-a))];
yv=zeros(length(y0),n+1); % matrix with row vectors
yv(:,1) = y0;

%SymplektiskEuler
for ii = 1:n
    r = sqrt(yv(1,ii)^2 + yv(2,ii)^2);
    yv(3,ii+1) = yv(3,ii) - h * (yv(1,ii) / r^3); %p1_n+1
    yv(4,ii+1) = yv(4,ii) - h * (yv(2,ii) / r^3); %p2_n+1

    yv(1,ii+1) = yv(1,ii) + h * yv(3,ii+1); %q1_n+1
    yv(2,ii+1) = yv(2,ii) + h * yv(4,ii+1); %q2_n+1
end


%plottar bana
figure;
plot(yv(1,:), yv(2,:), 'b');
xlabel('q1');
ylabel('q2');
title('Omloppsbana, symplektisk Euler');
grid on;
axis equal;

E = zeros(1,n+1); %energi 
for ii = 1:n+1
    r = sqrt(yv(1,ii)^2 + yv(2,ii)^2);
    E(ii) = 0.5 * (yv(3,ii)^2 + yv(4,ii)^2) - 1 / r;
end

figure;
plot(tv, E, 'r');
xlabel('Tid');
ylabel('Total energi');
title('Energins förändring över tid, symplektisk Euler');
grid on;
ylim([-0.51 -0.49])