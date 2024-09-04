function y = cmytail(ap,cp,mp,yp)
if yp>=0
    term1=(mp.*ap).^(-yp).*exp(-mp.*ap)./yp;
    term2=(mp.*ap).^(1-yp).*exp(-mp.*ap)./(yp.*(1-yp));
    term3=(gamma(2-yp)./(yp.*(1-yp))).*(1-gammainc(mp.*ap,2-yp));
    y=cp.*mp.^yp.*(term1+term2-term3);
else
    y=cp.*mp.^yp.*gamma(-yp).*(1-gammainc(mp.*ap,-yp));
end
return
