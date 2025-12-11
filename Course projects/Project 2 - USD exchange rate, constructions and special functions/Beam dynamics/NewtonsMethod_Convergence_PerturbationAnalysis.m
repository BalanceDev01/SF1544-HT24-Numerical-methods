%konstanter
%H = 4.1355432362; %För 2a) och 2b) ändra H=0.5 
H=0.5;
A = 9;
alpha = 1/2;
w = 4;
%Svängande balk, förskjutning från jämviktsläge 
y = @(t) A*exp(-alpha*t).*cos(w*t);

f=@(t) y(t) - H; %
fprim=@(t) -alpha*A*(exp(-alpha*t)).*cos(w*t) - A*exp(-alpha*t).*w.*sin(w*t);
%plottar funkt f jag vill hwitta nollställe till
subplot(2,1,1)
fplot(y,[-38 38]); 
hold on
plot([-38 38],[0 H], '--')
hold off
xlabel('x');
ylabel('f(x)');
legend('$y(x) = A  e^{-\alpha t} \cos(x) $', 'Interpreter', 'latex');
subplot(2,1,2) %behövs ändra plot för H = 4.1355432362 om man vill se rot
fplot(f,[4 6]); 
hold on
plot([4 6],[0 0], '--')
title('Sista nollstället tbar = 4.92731848') %för H=0.5 alltså
xlabel('x');
ylabel('f(x)');
legend('$y(x) = A  e^{-\alpha t} \cos(x) - 0.5$', 'Interpreter', 'latex');

%redan löst numeriskt och får "exakt" lösn t = 0.255374677
t_exakt = 4.92731848; %FÖR H=0.5 
%t_exakt = 1.53970693; %FÖR H=4.1355432362
% Parametrar för Newtons metod
diffx=1; iter=0; maxiter=100; %x är t

%tolerans
tol=10^(-8);
%startgissning
t_n = 5; %startgissn = 5 ger 4 iterationer och tydlig konv. för H=0.5
         %STARTGISSN = 5 FÖR H=4.1355432362 GER 21 ITER OCH ROT T_H = 1.539 

tvector=[t_n]; %vektor, sparar alla nollställen så jag kan ta sista för min givna tolerans
% Newtons metod ------------------------------------------------------
disp('t_(n+1)               e_n');
while diffx>tol&iter<maxiter
    iter = iter + 1;
    
    t_ny = t_n - f(t_n)/fprim(t_n);
   
    tvector(iter+1) = t_ny;
    
    % Beräkna noggranheten (x är t)
    diffx_ny=abs(tvector(iter+1)-tvector(iter)); %lite likt diffx = abs(t_ny - t_n);
    
    %p = log(diffx_ny)/log(diffx);%OPTIONAL:Beräknar ordningen av konvergens (kvadratisk)
    e_n= abs(t_ny - t_exakt);
    % Printa
    disp([t_ny,e_n]); 
    %disp(['konv.ordn. p = ', num2str(p)])%OPTIONAL
    % Updatera t_n
    t_n = t_ny;
    % Updaterar också diffx
    diffx = diffx_ny;%OPTIONAL
end
if (iter==maxiter) %Kollar ifall while-loop slutade pga max iterationer
    disp('========='); 
    disp('Warning: max no of iters reached.'); 
    disp('========='); 
end
if diffx <= tol %Vice versa,kollar denna ifall while-loop slutade pga konvergens
    disp('Lösningen konvergerade!');
else
    disp('Lösningen konvergerade ej');
end
%------------------------------------------------------------------------
disp(['Sista tid: t_bar = ', num2str(tvector(end),9),' , Funktionsvärde: f(t_bar) = ', num2str(f(tvector(end)))]);
disp(['Antal iterationer: ', num2str(iter)])

%kollar kravet för kvadratisk konvergens...
tbar=tvector(end);
ev=abs(tvector(1:end-1)-tbar);
disp('===='); 
disp('t_bar beskriver sista iteration.') 
disp('e_n=abs(t_n-tbar) för varje iteration:');
disp(ev(:));

disp('===='); 
disp('Klicka vadsomhelst för att se kvadratisk konvergens'); 
pause;
disp('===='); 
disp(''); 
disp('e(n+1)/e(n)^2:');
Mv=ev(2:end)./ev(1:end-1).^2;
disp(Mv(:))
disp(['Kan se att värde M = ', num2str(Mv(end)),' < oändlighet så kvadratisk konv.'])
                                                %FÖR H = 4.1355432362 går
                                                %typ faktiskt M --> OÄNDLIGHET

%%
%OBS, OM t = NaN GES SÅ KÖR KOD OVAN FÖRST. FÖRMODLIGEN EJ UPPDATERAT WORKSPACE
%OBS! GLÖM EJ BYTA H=0.5 FÖR 2b)
% parametrar
diffx = 1; iter = 0; maxiter = 100; 

%tolerans
tol = 10^(-8);
%startgissning
t_n = 5; t_n_0 = 6; %startgissn = 5 och = 6 ger bra antal iterationer och nollställe för H=0.5
disp(['Startgissn. t_n = ',num2str(t_n) ,' och t_n_0 = ',num2str(t_n_0) ,' .'])
tvector=[t_n];
% Sekantmetoden ------------------------------------------------------
%disp('t_(n+1)               abs(t_(n+1)-t_(n))'); %OPTIONAL
while diffx>tol&iter<maxiter
    iter = iter + 1;
    
    t_ny = t_n - f(t_n) * (t_n - t_n_0) / (f(t_n) - f(t_n_0));
   
    tvector(iter+1,1) = t_ny;

    %Beräkna noggranheten (x är t)
    diffx=abs(tvector(iter+1)-tvector(iter)); %lite likt diffx = abs(t_ny - t_n);
    
    % Printa 
    %disp([t_ny ,diffx]); %OPTIONAL

    %Updatera t_n och t_n_0
    t_n_0 = t_n;
    t_n = t_ny;
   
end
if (iter==maxiter) %Kollar ifall while-loop slutade pga max iterationer
    disp('========='); 
    disp('Warning: max no of iters reached.'); 
    disp('========='); 
end
if diffx <= tol %Vice versa,kollar denna ifall while-loop slutade pga konvergens
    disp('Lösningen konvergerade!');
else
    disp('Lösningen konvergerade ej');
end
%-------------------------------------------------------------
disp(['Sista tid: t_bar = ', num2str(t_n,9),' , Funktionsvärde: f(t_bar) = ', num2str(f(tvector(end)))]);
disp(['Antal iterationer: ', num2str(iter)])