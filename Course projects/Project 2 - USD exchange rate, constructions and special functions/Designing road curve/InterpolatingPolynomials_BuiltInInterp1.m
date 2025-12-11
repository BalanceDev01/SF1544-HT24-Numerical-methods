%Given data för punkter och längder
A_coord = [170 950; 420 2400; 670 1730]; % (x_A, y_A) för P1, P2, P3
B_coord = [160 1008; 370 2500; 640 1760]; % (x_B y_B) för P1, P2, P3
M = [A_coord, B_coord];
L_A_vect = [60; 75; 42]; % L_A för P1, P2, P3
L_B_vect = [45; 88; 57]; %L_B för P1, P2, P3
m = [L_A_vect, L_B_vect];

%finding a good starting guess
for i = 1:size(M, 1)
    %koordinater för denna iteration av punkt
    x_A = M(i, 1);
    y_A = M(i, 2);
    x_B = M(i, 3);
    y_B = M(i, 4);
    L_A = L_A_vect(i, 1);
    L_B = L_B_vect(i, 1);
    
    hold on;
    grid on;
    angle = linspace (0,2*pi,62.8);
    for j = linspace(1,3,2)
        % disp(m(i,1))
        % disp(m(i,2))
        % disp([num2str(M(i,j)),'och',num2str(M(i,j+1))])
        axis equal; 
        plot(m(i,1)*cos(angle)+M(i,j),m(i,1)*sin(angle)+M(i,j+1),'Color', 'b')
    end    
end
%GRAFISKT SES ATT:
start_gissn = [217 990; 437 2473; 681 1769];

%%
%redan räknat lösn av punkter så för att kolla e_n jämför jag "exakt" lösn
sol_exakt_list = [204.156642794528 999.328731516298; 346.552542575773 2415.18127128789; 696.957247632660 1762.20724763266];
%sparar lösn av punkter
sol_list = zeros(3,2);
% Newtons metod --------------------------------------
% itererar genom P1, P2, P3 genom att iterera genom raderna till M
for i = 1:size(M, 1)
    %koordinater för denna iteration av punkt
    x_A = M(i, 1);
    y_A = M(i, 2);
    x_B = M(i, 3);
    y_B = M(i, 4);
    L_A = L_A_vect(i, 1);
    L_B = L_B_vect(i, 1);
    %funkt. som måste definieras i for-loop för har olika konstanta värden för olika punkter 
    F1 = @(x_P, y_P, x_A, y_A, L_A) (x_P - x_A)^2 + (y_P - y_A)^2 - L_A^2;
    F2 = @(x_P, y_P, x_B, y_B, L_B) (x_P - x_B)^2 + (y_P - y_B)^2 - L_B^2;

    F = @(x_P, y_P) [F1(x_P, y_P, x_A, y_A, L_A); F2(x_P, y_P, x_B, y_B, L_B)]; 
    DF = @(x_P, y_P) [2*(x_P - x_A), 2*(y_P - y_A); 2*(x_P - x_B), 2*(y_P - y_B)];
    

    %parametrar
    tol = 1e-10; s = Inf; iter=0; maxiter=100;
    x_P = start_gissn(i,1);
    y_P = start_gissn(i,2);
    iter = 0;
    s = Inf;
        %Själva Newtonuppdateringen
        while norm(s) >= tol&iter<maxiter
            iter = iter + 1; %räknar vilken iteration vi är på
            %beräknar andra term s i taylorutveckling
            s = -DF(x_P, y_P) \ F(x_P, y_P);
            x_P = x_P + s(1);
            y_P = y_P + s(2);
            disp(['x: e_n = ', num2str(x_P - sol_exakt_list(i,1)) , ' , y: e_n = ',num2str(y_P - sol_exakt_list(i,2))])
            %disp(['Iteration: ',num2str(iter),'. x: e_n = ',num2str(norm(s(1))),' och y: e_n = ',num2str(norm(s(2))),'.']) % likt diffx då norm(s(1))= norm(x_P_n+1 - x_P_n)
        end
        %Kod nedan behövs ej! Eftersom grafiskt valt startgissn är så P
        %alltid fås, kod här startgissar så Q sen P ges och väljer P.
        %sol(k, :) = [x_P, y_P]; % 2x2-matris med lösn Q och P som kolumner
        sol_list(i,:) = [x_P, y_P];
    % %disp(sol) %Printar lösn P och Q till punkt för iteration 
    % if sol(1,1) > sol(2,1) %Kollar för just högra skärning P
    % sol_right = sol(1,:); 
    % else 
    % sol_right = sol(2,:);
    % end

    disp(['Lösn av koordinat till punkt P', num2str(i), ': (',num2str(sol_list(i,1)),',',num2str(sol_list(i,2)),') .'])
    
end

%----------------------------------------------------------------------------

%%
%POLYNOMISK MODELL
% i Ax=b med minstakvadratlösn c till A'Ac=A'b , y är  b
P1 = [0 0];
P5 = [1020 0];
x_vect = [P1(1,1); sol_list(:,1); P5(1,1)]; %x-värde på punkter P (OBS! BLANDA EJ IHOP MED x i workspace från ROADCOORD.MAT)
y_vect = [P1(1,2); sol_list(:,2); P5(1,2)]; %y-värde på punkter P (OBS! BLANDA EJ IHOP MED y i workspace från ROADCOORD.MAT)
xval = linspace(1,1020,1020); %BLANDA EJ IHOP MED xval FÖR LAGRANGEINTERPOL I 3c) 
% Vendermonde matris A:
A = [x_vect.^0, x_vect x_vect.^2 x_vect.^3 x_vect.^4];
c_poly = A\y_vect; %minstakvadratlösn
disp(['Konditionstalet är ',num2str(cond(A)) ,'.'])
disp('Värden på polynomkoeff. är: ')
disp(['c_0 = ',num2str(c_poly(1)) ,'.'])
disp(['c_1 = ',num2str(c_poly(2)) ,'.'])
disp(['c_2 = ',num2str(c_poly(3)) ,'.'])
disp(['c_3 = ',num2str(c_poly(4)) ,'.'])
disp(['c_4 = ',num2str(c_poly(5)) ,'.'])
%(Måste stoppa in xval nedan i funktion istället för x, vill ha kontinuerlig funktion inte bara punkter)
f_poly=@(xval,c_poly) c_poly(1) + c_poly(2).*xval + c_poly(3)*xval.^2 + c_poly(4)*xval.^3 + c_poly(5)*xval.^4 ; 
plot(xval,f_poly(xval,c_poly),'--','Color', 'r')
hold on;
plot(x_vect,y_vect,'o', 'MarkerSize', 10, 'Color', 'b')
xlabel('x_vect');
ylabel('y_vect');
grid on;
title('Bestämmning av vägsträckning');

%%
%(LAGRANGES INTERPOLATION*), (*) verkar ej ge helt rätt plot för 3d)! borde nog använda Newtons ansats istället
% i Ax=b med minstakvadratlösn c till A'Ac=A'b , y är  b
x = x; % gör inget, visar mest att x är från att ladda in ROADCOORD.MAT
y = y; % gör inget, visar mest att y är från att ladda in ROADCOORD.MAT 
%x = x2; % gör inget, visar mest att x2 är från ROADCOORD2.MAT  
%y = y2; % gör inget, visar mest att y2 är från ROADCOORD2.MAT 

deg=10; %ger en blå kurva som skär punkterna jag vill den ska skära
%deg=4; %samma grad som i 3b)

%plottar individuella punkter 
hold on;
plot(x,y,'o', 'MarkerSize', 5, 'Color', 'black')
xlabel('x');
ylabel('y');
grid on;
title('Annan vägsträckning');

xval = linspace(min(x), max(x));
for i=0:deg
    yval = lagrangepol(i, xval, x);
    plot(xval, yval)
end
grid on

%interpolerande polynom ges alltså av:
b = y;
polynomet = 0;
for i=0:deg
    polynomet = polynomet + b(i+1)*lagrangepol(i, xval, x);
end
plot(xval, polynomet, 'LineWidth', 2)
%Får RUNGES FENOMEN här i plotten med LAGRANGES INTERPOLATION! DETTA
%FÖRMODLIGEN ÄR PÅ GRUND AV ATT HÖGRE GRADENS POLYNOM SOM LAGRANGES
%INTERPOLATION ANVÄNDER GER STORA OSCILLATIONER. KAN FIXAS TILL GENOM
%ANVÄNDNING AV CHEBYSHEVSPUNKTER ISTÄLLET FÖR JÄMNT FÖRDELADE PUNKTER.
function varde = lagrangepol(i, xval, x) % xd = datapunkterna
    n = length(x) - 1; % datapunkterna heter x0, x1, ..., xn
    varde = 1;
    for k=0:n
        if i == k
            continue; % hoppar över denna iteration
        end
        varde = varde .* (xval - x(k+1)) ./ (x(i+1) - x(k+1)); % indexen börjar på 1, därför +1
    end
end



%%
%ANVÄNDER INTERP1

plot(x, y, 'o', 'MarkerSize', 5, 'DisplayName', 'Datapunkter');
hold on;

xx = linspace(-1,1,200);
v=interp1(x,y,xx);
plot(xx,v,'r')

%DÅ INTERP1 SOM DEFAULT ANVÄNDER LINJÄR INTERPOLATION (KOPPLAR BARA IHOP
%PUNKTER) SÅ ÄR POLYNOMGRAD HÄR EJ HÖG OCH RUNGES FENOMEN UPPSTÅR EJ.

%%
% NEWTONS ANSATS (Använd clear och ladda sen in variabler på nytt)
x = x'; %transponerar bara, men visar att x laddas in från roadcoord.mat
y = y'; %transponerar bara, men visar att y laddas in från roadcoord.mat
x2 = x2'; %transponerar bara, men visar att x2 laddas in från roadcoord2.mat
y2 = y2'; %transponerar bara, men visar att y2 laddas in från roadcoord2.mat

clf;
plot(x, y, 'o', 'MarkerSize', 5, 'DisplayName', 'Datapunkter','Color','b');
hold on;
plot(x2, y2, 'o', 'MarkerSize', 5, 'DisplayName', 'Datapunkter','Color','g');
hold on;

Q = [x,y,x2,y2];
for ii = linspace (1,3,2)
    %räkn. Newtons delade differenser
    n = length(Q(:,ii)); %Q(:,ii) = x för ii = 1 sen Q(:,ii) = x2 för ii = 3
    div_diff = zeros(n, n);
    div_diff(:, 1) = Q(:,ii+1);%Q(:,ii+1) = y för ii = 1 sen Q(:,ii+1) = y2 för ii = 3
    
    %fyller tabell med delade diff
    for j = 2:n
        for i = j:n
            div_diff(i, j) = (div_diff(i, j-1) - div_diff(i-1, j-1)) / (Q(i,ii) - Q(i-j+1,ii)); 
        end                                        % Q(i,ii) - Q(i-j+1,ii) = x(i) - x(i-j+1)
    end
    
    %ta koeff för Newtons interpol.polynom
    coeffs = diag(div_diff);
    
    %checkar interpolationspolynomet vid punkter
    xval = linspace(min(Q(:,ii)), max(Q(:,ii)), 100);
                  % min(x)          max(x)
    P = coeffs(1) * ones(size(xval)); %kör först konstantledet
    
    for k = 2:n
        P = P + coeffs(k) .* prod(xval - Q(1:k-1,ii), 1);
    end                                % x(1:k-1)
    plot(xval, P, 'LineWidth', 2);
end
%plots
grid on;
xlabel('x');
ylabel('y');
title('Newtons Interpolationspolynomer');
legend show;

clear; %rensar variabler för annars fick jag problem med dimensionerna andra gången jag körde program

%Får RUNGES FENOMEN här i plotten med NEWTONS ANSATS för roadcoord.mat. 
% DETTA FÖRMODLIGEN ÄR PÅ GRUND AV ATT HÖGRE GRADENS POLYNOM SOM NEWTONS ANSATS
% ANVÄNDER GER STORA OSCILLATIONER VID ÄNDAR, ORIMLIG VÄG.

% MEN, FÖR roadcoord2.mat ÄR PUNKTER EJ HELT JÄMNT FÖRDELADE UTAN 
% FÖRDELNINGEN VERKAR STÄMMA ÖVERENS MED CHEBYSHEVPUNKTER SOM T.EX. GÖR
% ATT VID ÄNDARNA STÄLLS FLERA PUNKTER IN SÅ ATT OSCILLATIONERNA EJ UPPSTÅR