function [tv, yv] = implicit_trapezoid(A,g,tspan,y0,n)
    h=(tspan(2)-tspan(1))/n;
    tv=(tspan(1)+h*(0:n));
    p=length(y0);
    yv=zeros(p,n+1);
    yv(:,1)=y0;

    I = eye(size(A));
    M = I - (h/2)*A; % Implicita matrisen
    
    for ii=1:n
        g_n = g(tv(ii));
        g_np1 = g(tv(ii+1));
        rhs = yv(:,ii) + (h/2)*(A*yv(:,ii) + g_n + g_np1);
        yv(:,ii+1) = M\rhs;
    end
end
