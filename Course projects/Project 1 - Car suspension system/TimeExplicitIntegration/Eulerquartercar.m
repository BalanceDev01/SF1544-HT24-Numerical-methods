function [tv,yv]=Eulerquartercar(f,tspan,y0,n)
%Solve y'=f(t,y) with inital condition y(tspan(1))=y0
%up to time tspan(2).
%y0 should be a column vector,
%f is assumed to return a column vector of the same length.

h=(tspan(2)-tspan(1))/n;
tv=(tspan(1)+h*(0:n)); %tv a row vector.
p=length(y0);
yv=zeros(p,n+1);
yv(:,1)=y0;

for ii=1:n
    yv(:,ii+1)=yv(:,ii)+h*f(tv(ii),yv(:,ii));
end
