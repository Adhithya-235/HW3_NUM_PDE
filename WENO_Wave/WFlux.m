function [FL] = WFlux(xv,uv,h,beta,Crec,dw,m)

nx = length(xv);
[xe,ue] = extend(xv,uv,h,m);
um = zeros(nx+2,1);
up = zeros(nx+2,1);
for i = 1:nx+2
   
    [um(i),up(i)] = WENO(xe(i:(i+2*(m-1))),ue(i:(i+2*(m-1))),m,Crec,dw,beta);
    
end

FL = (-1/h)*(LaxF(up(2:nx+1),um(3:nx+2),1) - LaxF(up(1:nx),um(2:nx+1),1));

end