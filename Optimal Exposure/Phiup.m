function phi = Phiup(a,c,gam,lam)
    L = (lam-1)*(1+gam)/(a*c);
    f = @(u)( ( u / ( (1-u).^(gam/(1+gam)) ) - L )^2 );
    if lam < 1
        phi = 0;
    elseif lam == 1
        phi = a;
    elseif lam>1
        options = optimoptions('fmincon','Display','Off');
        u = fmincon(f,0.5,[],[],[],[],0,1,[],options);
        phi = -(1-lam).*log(u)/c+a*(1-u).^(1/(1+gam));
        %fprintf('u = %d, f(u) = %d\n', u,f(u))
    end
end