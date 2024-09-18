function p = InverseContract(cp,M,yp,cn,G,yn,xp,xn,x,delta)
    nu = [cn*exp(G*xn).*(-xn).^(-1-yn),cp*exp(-M*xp).*xp.^(-1-yp)]*delta;
    p = nu*(exp(-x')-1);
end