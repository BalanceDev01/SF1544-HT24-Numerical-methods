%Trapetsregeln för räkn av integral
color = ['b' 'r' 'k'];
xlist = [0.11 0.32 1.14];
for i = 1:3
    f = @(t) (2/sqrt(pi)) .* exp(-t.^2); % integranden
    f_bis = @(t) (4/sqrt(pi)).*abs(exp(-t^2)*(2*t^2-1));
    b= xlist(1,i); %gör sen loop för alla 3 x-värden
    a = 0;
    % (b) Beräkning av värdet
    N_range = 20:700;
    N = 100 ; % antal delintervall alt: räkna N = (b-a)/h med vald h, BYT DÅ VILKEN RAD KODEN ÄR PÅ 
    h = (b-a)/N ; % alt: räkna h=(b-a)/N med vald N
    xv = a:h:b; % x-värden
    yv = f(xv); % y-värden
    Q = h*(yv(1)/2 + yv(end)/2 + sum(yv(2:end-1)))
    
    

    %felgräns 
    C = @(x) (1/12).*abs(f_bis(x));
    %ONÖDIGT, FUNKAR EJ HELT
    % if xlist(1,i) <= 1/sqrt(2) %maximum ges då x = 1/sqrt(2) så innan det ges f_max av f(x), mellan 0 till x
    %     C = @(x) (1/12).*abs(f_bis(x));
    % elseif xlist(1,i) > 1/sqrt(2) 
    %     C = @(x) (1/12).*abs(f_bis(1/sqrt(2))); 
    % end
    E_N_max =@(N_range, x) C(x).*(x^3)./(N_range.^2);
    plot(N_range,E_N_max(N_range,xlist(1,i)), "Color",color(i),'DisplayName', sprintf('E_N [g], x= %.2f', xlist(1,i)))
    hold on;
    
end
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlim([10^1, 10^3]);
ylim([10^-10, 10^-3])
grid on;
hold off;
legend show;
%-------------------------------------------
%%
f = @(t) (2/sqrt(pi)) * exp(-t.^2); 
color = ['b' 'r' 'k']; 
Nlist = [50 120 400]; 

for i = 1:3
    E_N = [];
    xval = linspace(0, 6, Nlist(i)); 
    
    for ii = 1:length(xval)
        b = xval(ii);
        a = 0;
        N = Nlist(i); 
        h = (b-a)/N; 
        
        
        xv = linspace(a, b, N+1); 
        yv = f(xv);
        Q = h * (0.5 * yv(1) + sum(yv(2:end-1)) + 0.5 * yv(end)); %approx. av integral
        
        %felet
        E_N(ii) = abs(Q - erf(b));
    end
    
    
    plot(xval, E_N, "Color", color(i), 'DisplayName', sprintf('E_N [g], N= %.0f', Nlist(i)));
    hold on;
end
xlim([0, 6.5]);
ylim([10^-16, 10^-4])
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;
legend show;
hold off;
