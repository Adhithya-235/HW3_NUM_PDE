clear
close all
clc

%PROBLEM DEFINITION
xs = -1;
xe = 1;
h = 0.020;
te = 2;
nx = 2/h;
x = (xs:h:xe)';
x_aux = (xs:0.5*h:xe)';
x_fv = x_aux(2:2:201);
W = zeros(length(x),1);
for i = 1:length(x)
   W(i) = wavetest(x(i));
end
%u0 = W(2:2:end);
u0 = W;
m = 5; %For a 2m-1 order scheme

%COMPUTE RECONSTRUCTION WEIGHTS
Crec = zeros(m+1,m);
for r = -1:m-1
   
    Crec(r+2,:) = ReconstructWeights(m,r);
    
end

%LINEAR COMBINATION WEIGHTS 
dw = LinearWeights(m,0);

%SMOOTHNESS INDICATOR
beta = zeros(m,m,m);
for r = 0:m-1
   
    xl = -0.5 + (-r:1:m-r);
    beta(:,:,r+1) = betarcalc(xl,m);
    
end

%TIMESTEPPING
it = 0;
t = 0;
u_old = u0;
CFL = 0.2;

while(t<te)
   
    
    aa = max(abs(u_old));
    k = CFL*h/aa; %ADAPTIVE TIMESTEPPING
    if ((t+k)>te)
       k = te-t; 
    end
    
    FL = WFlux(x,u_old,h,beta,Crec,dw,m);
    u = u_old + k*FL;
    FL = WFlux(x,u,h,beta,Crec,dw,m);
    u = 0.25*(3*u_old + u + k*FL);
    FL = WFlux(x,u,h,beta,Crec,dw,m);
    u = (1/3)*(u_old + 2*u + 2*k*FL);
    
    u_old = u;
    plot(x,u,'r',x,u0,'k','linewidth',2)
    axis([-1 1 -0.2 1.2])
    drawnow
    
    it = it+1;
    t = t+k;
    
end


