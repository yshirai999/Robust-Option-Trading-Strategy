function c = NA(cp,M,etap,yp,cn,G,etan,yn,xp,xn,x,delta)
    p = InverseContract(cp,M,yp,cn,G,yn,xp,xn,x,delta);
    M = M+etap;
    G = G+etan;
    nu = [cn*exp(G*xn).*(-xn).^(1-yn),cp*exp(-M*xp).*xp.^(1-yp)]*delta;    
    c = (nu*(exp(-x')-p)).^2+(nu*(exp(x')-1)).^2;
end