function phi = Phitil(b,c,lam)
    lamp = lam(lam<1);
    phi(lam<1) = -(1-lamp).*(log((1-lamp)/b)/c)-(b/c)*(1 - (1-lamp)/b );
    phi(lam==1) = -b/c;
    phi(lam==0) = (log(b)+1-b)/c;
end