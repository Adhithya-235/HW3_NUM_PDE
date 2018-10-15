function[um,up] = WENO(xloc,uloc,m,Crec,dw,beta)

p = 1; 
q = m-1; 
vareps = 1e-6;

if m==1
    um = uloc(1); up = uloc(1);
else
    alpham = zeros(m,1); 
    alphap = zeros(m,1);
    upl = zeros(m,1); 
    uml = zeros(m,1);
    betar = zeros(m,1);
end

for r = 0:m-1
    umh = uloc(m-r+[0:m-1]);
    upl(r+1) = Crec(r+2,:)*umh;
    uml(r+1) = Crec(r+1,:)*umh;
    betar(r+1) = (umh')*beta(:,:,r+1)*umh;
end

alphap = dw./((vareps + betar).^(2*p));
alpham = flipud(dw)./((vareps + betar).^(2*p));

um = (alpham')*uml/sum(alpham);
up = (alphap')*upl/sum(alphap);

end


