clear 
close all
clc

%PROBLEM DEFINITION
veps = 1e-6;
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
u = u0;
up = zeros(size(u));

%RECONSTRUCTION FUNCTIONS

    
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
        
        
            
        
    end
end