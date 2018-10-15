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
x_grid = (xs:dx:xe)';
x_aux = (xs:0.5*dx:xe)';
x_fv = x_aux(2:2:end);
W = 0.25 + 0.5*sin(pi*x_grid);
%u0 = W(2:2:end);
u0 = W;
TV = zeros(nt,1);

for m = 1:4
    u = u0;
    up = zeros(size(u));
    for it = 1:nt
        alph = max(u);
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
            
            flfp = 0.25*(u(k)^2+u(kp)^2) - 0.5*alph*(u(kp)-u(k));
            flfm = 0.25*(u(k)^2+u(km)^2) - 0.5*alph*(u(k)-u(km));
            flwp = 0.25*(u(k)^2+u(kp)^2) - 0.125*lam*(u(kp)+u(k))*(u(kp)^2-u(k)^2);
            flwm = 0.25*(u(k)^2+u(km)^2) - 0.125*lam*(u(km)+u(k))*(u(k)^2-u(km)^2);
            if u(k)>=0 
                rp = (u(k) - u(km))/(u(kp) - u(k));
                rm = (u(km) - u(kmm))/(u(k) - u(km));
            elseif u(k)<=0
                rp = (u(kpp) - u(kp))/(u(kp) - u(k));
                rm = (u(kp) - u(k))/(u(k) - u(km));
            end
            fbwp = flfp + rp*(flwp-flfp);
            fbwm = flfm + rm*(flwm-flfm);
            ffp = 0.5*(fbwp + flwp);
            ffm = 0.5*(fbwm + flwm);
            
            switch m
                case 1 %UPWIND
                    up(k) = u(k) - lam*(flfp - flfm);
                case 2 %LAX WENDROFF
                    up(k) = u(k) - lam*(flwp-flwm);
                case 3 %BEAM WARMING
                    up(k) = u(k) - lam*(fbwp-fbwm);
                case 4 %FROMM's METHOD
                    up(k) = u(k) - lam*(ffp-ffm);                    
            end                            
        end
        u = up;
        figure(1)
        subplot(2,2,m)
        plot(x_grid,u,'r',x_grid,u0,'k','linewidth',2)
        axis([-1 0 -0.4 0.8])
        drawnow
        
        TV(it) = sum(abs(u(2:nx)-u(1:nx-1)));
    end
    
    figure(2)
    subplot(2,2,m)
    loglog(1:nt,TV,'k','linewidth',2)
    title('Trend in Total Variation')
    xlabel('Timestep')
    ylabel('TV')
    %axis([-1 0 -0.4 0.8])
    %drawnow
        

end
