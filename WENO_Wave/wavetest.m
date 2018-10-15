function[W] = wavetest(x)

a = 0.5;
z = -0.7;
dta = 0.005;
alp = 10;
bet = log(2)/(36*dta*dta);

G = @(x,b,ze) exp(-b*((x - ze)^2));
F = @(x,al,aa) sqrt(max((1-al*al*((x-aa).^2)),0));

if ((x>=-0.8) && (x<=-0.6))
    W = (1/6)*(G(x,bet,z-dta) + G(x,bet,z+dta) + 4*G(x,bet,z));
elseif ((-0.4<=x) && (x<=-0.2))
    W = 1;
elseif ((0<=x) && (x<=0.2))
    W = 1-abs(10*(x-0.1));
elseif ((0.4<=x) && (x<=0.6))
    W = (1/6)*(F(x,alp,a-dta) + F(x,alp,a+dta) + 4*F(x,alp,a));
else
    W = 0;
end

    
end