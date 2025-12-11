%Given data för punkter och längder
A_coord = [170 950; 420 2400; 670 1730]; % (x_A, y_A) för P1, P2, P3
B_coord = [160 1008; 370 2500; 640 1760]; % (x_A, y_A) för P1, P2, P3
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
    xlabel('x')
    ylabel('y')
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
    coord_vect=[x_P,y_P];
    iter = 0;
    s = Inf;
        %Själva Newtonuppdateringen
        while norm(s) >= tol&iter<maxiter
            iter = iter + 1; %räknar vilken iteration vi är på
            %beräknar andra term s i taylorutveckling
            s = -DF(x_P, y_P) \ F(x_P, y_P);
            x_P_ny = x_P + s(1);
            y_P_ny = y_P + s(2);
            coord_vect(iter+1,1) = x_P_ny;
            coord_vect(iter+1,2) = y_P_ny;
            %uppdaterar
            x_P = x_P_ny;
            y_P = y_P_ny;
            %disp(['x: e_n = ', num2str(x_P_ny - sol_exakt_list(i,1)) , ' , y: e_n = ',num2str(y_P_ny - sol_exakt_list(i,2))])
            %disp(['Iteration: ',num2str(iter),'. x: e_n = ',num2str(norm(s(1))),' och y: e_n = ',num2str(norm(s(2))),'.']) % likt diffx då norm(s(1))= norm(x_P_n+1 - x_P_n)
        end
        %Kod nedan behövs ej! Eftersom grafiskt valt startgissn är så P
        %alltid fås, kod här startgissar så Q sen P ges och väljer P.
        %sol(k, :) = [x_P, y_P]; % 2x2-matris med lösn Q och P som kolumner
        sol_list(i,:) = [x_P_ny, y_P_ny];
    % %disp(sol) %Printar lösn P och Q till punkt för iteration 
    % if sol(1,1) > sol(2,1) %Kollar för just högra skärning P
    % sol_right = sol(1,:); 
    % else 
    % sol_right = sol(2,:);
    % end

    disp(['Lösn av koordinat till punkt P', num2str(i), ': (',num2str(sol_list(i,1)),',',num2str(sol_list(i,2)),') .'])
    
end
%----------------------------------------------------------------------------
%kollar kravet för kvadratisk konvergens...
xbar=coord_vect(end,1);
ybar=coord_vect(end,2);
ev_x=abs(coord_vect(1:end-1,1)-xbar);
ev_y=abs(coord_vect(1:end-1,2)-ybar);
disp('===='); 
disp('x_bar och y_bar beskriver sista iteration i x och y.') 
disp('e_n=abs(x_n-xbar) och e_n=abs(y_n-ybar) för varje iteration:');
disp([ev_x(:),ev_y(:)]);

disp('===='); 
disp('Klicka vadsomhelst för att se kvadratisk konvergens'); 
pause;
disp('===='); 
disp(''); 
disp('x: e(n+1)/e(n)^2 och y: e(n+1)/e(n)^2');
Mv_x=ev_x(2:end,1)./ev_x(1:end-1,1).^2;
Mv_y=ev_y(2:end,1)./ev_y(1:end-1,1).^2;
disp([Mv_x(:), Mv_y(:)])
disp(['Kan se att värde M_x = ', num2str(Mv_x(end-1)),' och M_y = ',num2str(Mv_y(end-1)),' < oändlighet så kvadratisk konv.'])

%%
%POLYNOMISK MODELL
% i Ax=b med minstakvadratlösn c till A'Ac=A'b , y är  b
P1 = [0 0];
P5 = [1020 0];
x = [P1(1,1); sol_list(:,1); P5(1,1)]; %x-värde på punkter P
y = [P1(1,2); sol_list(:,2); P5(1,2)]; %y-värde på punkter P
xval = linspace(1,1020,1020);
% Vendermonde matris A:
A = [x.^0, x x.^2 x.^3 x.^4];
c_poly = A\y; %minstakvadratlösn
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
plot(x,y,'o', 'MarkerSize', 10, 'Color', 'b')
xlabel('x');
ylabel('y');
grid on;
title('Bestämmning av vägsträckning');