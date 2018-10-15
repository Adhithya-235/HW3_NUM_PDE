function [xe,ue] = extend(x,u,h,m)
% Purpose: Extend dependent and independent vectors (x,u) , by m cells % subject to approproate boundary conditions .
% BC = "D" ? Dirichlet; BC = "N" ? Neumann; BC = "P" ? periodic
% ul/ur ? BC value ? only active for Dirichlet BC
xl=min(x); xr=max(x);N=length(u);
xe = zeros(N+2*m,1); ue = zeros(N+2*m,1); q = 1:m;
% Extend x
xe(m-q+1) = xl-q*h; xe(N+m+q) = xr + q*h; 
xe((m+1):(N+m)) = x(1:N);
ue(m-q+1) = u(N-q); 
ue(N+m+q) = u(q+1); 
ue((m+1):(N+m)) = u(1:N); 
end
