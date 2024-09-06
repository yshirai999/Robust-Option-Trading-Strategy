a=importdata('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\spyccgmyy2sommats20231231output.mat');
days=unique(a(:,1));
ndays=length(days);
parm=abs(a(:,5:10));
cpv=parm(:,1);
cnv=parm(:,2);
bpv=1./parm(:,4);
bnv=1./parm(:,3);
ypv=2*exp(-abs(parm(:,5)));
ynv=2*exp(-abs(parm(:,6)));
parmm=[cpv bpv ypv cnv bnv ynv];
























