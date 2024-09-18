function [vv,ss] = NUG(pv,uv,bv,var,N)
s=log(pv);
b=log(uv);
a=log(bv);
gg=var;
g1=gg;
g2=(b-s)/(s-a)*g1;
c1=asinh((a-s)./g1);
c2=asinh((b-s)./g2);
xv=zeros(N/2,1);
yv=zeros(N/2,1);
for ik=1:N/2
    x=s+g1.*sinh(c1.*(1-(ik-1)./(N/2-1)));
    y=s+g2.*sinh(c2.*2.*ik./N);
    xv(ik)=x;
    yv(ik)=y;
end
vv=[xv' yv']';
ss=exp(vv);
return