clear 
close all
clc

%PROBLEM DEFINITION
xs = -1;
xe = 1;
dx = 0.020;
te = 2;
dt = 0.018;
nx = 2/dx;
nt = ceil(te/dt);
lam = (dt/dx);
x_aux = (xs:0.5*dx:xe)';
x_fv = x_aux(2:2:201);
W = zeros(length(x_aux),1);
for i = 1:length(x_aux)
   W(i) = wavetest(x_aux(i));
end
u0 = W(2:2:end);
ErrMat = zeros(4,2);

for m = 1:4
    u = u0;
    up = zeros(size(u));
    for it = 1:nt
        for k = 1:nx
            
            kp = k + 1;
            kpp = k + 2;
            km = k - 1;
            kmm = k - 2;
            
            if (kp > nx)
                kp = kp - (nx);
            end
            if (kpp > nx)
                kpp = kpp - (nx);
            end
            if (km < 1)
                km = km + (nx);
            end
            if (kmm < 1)
                kmm = kmm + (nx);
            end
            
            flfp = u(k);
            flfm = u(km);
            flwp = 0.5*(u(k)+u(kp)) - 0.5*lam*(u(kp)-u(k));
            flwm = 0.5*(u(k)+u(km)) - 0.5*lam*(u(k)-u(km));
            rp = (u(k) - u(km))/(u(kp) - u(k));
            rm = (u(km) - u(kmm))/(u(k) - u(km));
            
            switch m
                
                case 1 %MINMOD
                    fp = flfp + (max([0,min(rp,1)]))*(flwp - flfp);
                    fm = flfm + (max([0,min(rm,1)]))*(flwm - flfm);
                    up(k) = u(k) - lam*(fp-fm);
                case 2 %SUPERBEE
                    fp = flfp + (max([0,min(2*rp,1),min(rp,2)]))*(flwp - flfp);
                    fm = flfm + (max([0,min(2*rm,1),min(rm,2)]))*(flwm - flfm);
                    up(k) = u(k) - lam*(fp-fm);
                case 3 %MC
                    fp = flfp + (max([0,min([2*rp,0.5*(1+rp),2])]))*(flwp - flfp);
                    fm = flfm + (max([0,min([2*rm,0.5*(1+rm),2])]))*(flwm - flfm);
                    up(k) = u(k) - lam*(fp-fm);
                case 4 %VAN LEER
                    fp = flfp + (max([0,(2*rp)/(1+rp)]))*(flwp - flfp);
                    fm = flfm + (max([0,(2*rm)/(1+rm)]))*(flwm - flfm);
                    up(k) = u(k) - lam*(fp-fm);                    
            end
        end
        u = up;
        figure(1)
        subplot(2,2,m)
        plot(x_fv,u,'r',x_fv,u0,'k','linewidth',2)
        axis([-1 1 -0.2 1.2])
        drawnow
    end
    ErrMat(m,1) = norm(u-u0,1);
    ErrMat(m,2) = norm(u-u0,2);
end
