% Datapunkter för dollarkurs under N dagar (730dagar=2år):
N=730;
L=400;
days = linspace(1,N,N)';
dvalue = USDSEK;

%LINJÄR MODELL
% i Ax=b med minstakvadratlösn c till A'Ac=A'b är dvalue b
% linjär matris A:
A_linj = [days.^0, days];
c_linj = A_linj\dvalue; %minstakvadratlösn
c_linj_new = (A_linj'*A_linj)\(A_linj'*dvalue); %dålig minstakvadratlösn
disp(['Värdena på linj. koeff. är c0 = ',num2str(c_linj(1)), ', c1 = ' ,num2str(c_linj(2)), '.'])
disp(['Skillnad i koefficienter Δc0 = ',num2str(c_linj(1)-c_linj_new(1)), ' och Δc1 = ' ,num2str(c_linj(2)-c_linj_new(2)), '.']) %skillnad mellan dålig och bra minstakvadraträkn
disp(['Konditionstalet är ',num2str(cond(A_linj)) ,'.'])
f_linj=@(t,y) (c_linj(1) + c_linj(2)*t);
MSE_linj = (1/N)*(norm(dvalue - f_linj(days,dvalue)).^2); %Samma som MSE = (1/N)*sum((dvalue - f(days,dvalue)).^2)
errv_linj = abs(dvalue - f_linj(days,dvalue));

%TRIGONOMETRISK MODELL
% i Ax=b med minstakvadratlösn c till A'Ac=A'b är dvalue b
% Vendermonde matris A:
A = [days.^0, days sin(((2*pi)/L)*days) cos(((2*pi)/L)*days)];
d = A\dvalue; %minstakvadratlösn
% d_new = (A'*A)\(A'*dvalue); %långsammare sätt att räkna d 
disp('Värdena på trig. koeff. är: ')
disp(['d_0 = ',num2str(d(1)) ,'.'])
disp(['d_1 = ',num2str(d(2)) ,'.'])
disp(['d_2 = ',num2str(d(3)) ,'.'])
disp(['d_3 = ',num2str(d(4)) ,'.'])
disp(['L = ',num2str(L) ,'.'])
f=@(t,y) d(1) + d(2)*t + d(3)*sin((2*pi/L)*t) + d(4)*cos((2*pi/L)*t);

MSE = (1/N)*(norm(dvalue - f(days,dvalue)).^2); %Samma som MSE = (1/N)*sum((dvalue - f(days,dvalue)).^2)

errv = abs(dvalue - f(days,dvalue)); %felet, (mätdata vs anpassad funktion)

%GAUSS-NEWTONS METOD i flera variabler R^n till R^m -------------------
c = [d(1);d(2);d(3);d(4);L ]; % startgissning
f_ny=@(t,y,c) c(1) + c(2)*t + c(3).*sin((2.*pi/c(5)).*t) + c(4).*cos((2.*pi/c(5)).*t);
r= @(t,y,c) f_ny(t,y,c)-y; %residualen
D_r = @(t, y, c) [t.^0, t, sin((2.*pi/c(5)).*t), cos((2.*pi/c(5)).*t), -c(3).*(2.*pi.*t/c(5)^2).*cos((2.*pi/c(5)).*t)+c(4).*(2.*pi.*t/c(5)^2).*sin((2.*pi/c(5)).*t)];
tol = 1e-10; s = Inf; iter=0; maxiter=100;
while norm(s) >= tol&&iter<maxiter
    s = -D_r(days,dvalue,c) \ r(days,dvalue,c); % överbestämt system!
    %disp(norm(s)) % kan användas för att studera metodens konvergens
    c = c + s;
end
if (iter==maxiter) %Kollar ifall while-loop slutade pga max iterationer
    disp('Warning: max no of iters reached.');
end
if norm(s) <= tol %Vice versa,kollar denna ifall while-loop slutade pga konvergens
    %disp('Lösningen konvergerade!');
else
    %disp('Lösningen konvergerade ej');
end
disp('Optimala värden (GN) på koeff.: ')
disp(['d0= ',num2str(c(1)), '.'])
disp(['d1= ',num2str(c(2)), '.'])
disp(['d2= ',num2str(c(3)), '.'])
disp(['d3= ',num2str(c(4)), '.'])
disp(['L= ',num2str(c(5)), '.'])
errv_ny = abs(dvalue - f_ny(days,dvalue,c));
%-----------------------------------------------------------------------


MSE_ny = (1/N)*(norm(dvalue - f_ny(days,dvalue,c)).^2); %Samma som MSE = (1/N)*sum((dvalue - f(days,dvalue)).^2)
disp(['Medelkvadratfelet för linjär anpassning är ', num2str(MSE_linj),'.'])
disp(['Medelkvadratfelet för trigonometrisk anpassning är ', num2str(MSE),'.'])
disp(['Medelkvadratfelet för Gauss-Newton är ', num2str(MSE_ny),'.'])
disp('  ')
disp(['Enligt linj. anpassning ändrades dollarkurs genomsnittliga värde med ',num2str(c_linj(2)),' per dag.'])
disp(['Enligt trig.anpassning ändrades dollarkurs genomsnittliga värde med',num2str((f(days(end),dvalue(end))-f(days(1),dvalue(1)))/N),' per dag.'])
disp(['Enligt GN-anpassning ändrades ddollarkurs genomsnittliga värde med',num2str((f_ny(days(end),dvalue(end),c)-f_ny(days(1),dvalue(1),c))/N),' per dag.'])
%plots
plot (days, f_linj(days,dvalue)) 
hold on;
plot (days, f(days,dvalue)) 
hold on;
plot (days, f_ny(days,dvalue,c)) 
hold on;
plot(days, dvalue, '.', color='g') %Det här är mätdatan
legend("Linjär anpassning","Trig. anpassning", "Gauss-Newt. anpassning");
xlabel('Tid[dagar]')
ylabel('Värde [SEK/USD]')
title('Anpassning')