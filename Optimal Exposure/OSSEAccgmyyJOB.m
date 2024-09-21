clear
clc
close all

%% Data

a=importdata('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\Data\spyccgmyy2sommats20231231output.mat');
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

alpha=1.2;
data = load('DSPEAMinputsnug.mat');
DSPEAMinputs1nug=data.DSPEAMinputs1nug;
DSPEAMinputs2nug=data.DSPEAMinputs2nug;
% for id=1:ndays
%     disp(id);
%     par1=parmm(1+2*(id-1),:);
%     par2=parmm(2+2*(id-1),:);
%     [matp1,matn1]=OSSEAccgmyyMeasures(par1,alpha);
%     [matp2, matn2]=OSSEAccgmyyMeasures(par2,alpha);
%     DSPEAMinputs1nug(id).measuresp=matp1;
%     DSPEAMinputs1nug(id).measuresn=matn1;
%     DSPEAMinputs2nug(id).measuresp=matp2;
%     DSPEAMinputs2nug(id).measuresn=matn2;
% end
% save('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\Optimal Exposure\DSPEAMinputsnug.mat','DSPEAMinputs1nug','DSPEAMinputs2nug');

%% Distortion
lamp=(1:.002:2);
lamn=(0:.002:1.0);
nlp=length(lamp); 
nln=length(lamn);

a=2; b=1; c=.5; gm=.75;
Phip = zeros(1,nlp);
for i=1:nlp
    Phip(i) = Phiup(a,c,gm,lamp(i));
end
Phin = -Phitil(b,c,lamn);
% Phip=zeros(1,nlp);
% for ic=1:nlp
%       lda=lamp(ic);
%       val=MeasureDistortionGPlusDualPhi(lda,[b c gm]);
%       Phip(ic)=val;
% end
% Phin=MeasureDistortionGMinusDualPhi(lamn,[b c gm]);
% Phin(nln)=b/c;

%% Exposure maximization
xx=nonzeros((-.5:.02:.5));
nx=length(xx);
verbose='False';
ExposureAcceptability=struct;
ndays_job = 10;
for id=1:ndays_job
      %disp(id);
      matp1=DSPEAMinputs1nug(id).measuresp;
      matn1=DSPEAMinputs1nug(id).measuresn;
      matp2=DSPEAMinputs2nug(id).measuresp;
      matn2=DSPEAMinputs2nug(id).measuresn;
      xxxxp=matp1(:,1);
      xxxxn=matn1(:,1);
      yyap=matp1(:,2);
      yyan=matn1(:,2);
      yy1p=matp1(:,3);
      yy1n=matn1(:,3);
      yy2p=matp1(:,4);
      yy2n=matn1(:,4);
      yy0p=matp1(:,5);
      yy0n=matn1(:,5);
      qqp=matp2(:,3);
      qqn=matn2(:,3);
      nxp=length(xxxxp);
      nxn=length(xxxxn);
      nxx=nxp+nxn;
      % yy=xx;
      xxxxx=[xxxxn' xxxxp']';
      A=zeros(nxx,nx);
      ee=eye(nx);
      for ic=1:nx
            A(:,ic)=interp1(xx,ee(:,ic),xxxxx,'spline','extrap');
      end
      % zz=.1*randn(1,nxx);
           
      % yyy=A*yy;
      
      wwa=[yyan' yyap']';
      ww1=[yy1n' yy1p']';
      ww2=[yy2n' yy2p']';
      ww0=[yy0n' yy0p']';
      wwq=[qqn' qqp']';

      xxxxx2=xxxxx.^2;

      ky=nx;
      kz=nxx;
   
      xxxxx2b=1./xxxxx2;
      pyenv('Version',...
          'C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\maximin2\python.exe',...
          'ExecutionMode','OutOfProcess');
      py.importlib.import_module('numpy');

      wwwa = py.numpy.array(wwa.');
      www1 = py.numpy.array(ww1.');
      www2 = py.numpy.array(ww2.');
      www0 = py.numpy.array(ww0.');
      wwwq = py.numpy.array(wwq.');
      lampp = py.numpy.array(lamp.');
      lamnn = py.numpy.array(lamn.');
      Phipp = py.numpy.array(Phip.');
      Phinn = py.numpy.array(Phin.');
      AA = py.numpy.array(A.');
      xxxxxa = py.numpy.array(xxxxx.');
      xxxxx2a= py.numpy.array(xxxxx2.');
      xxxxx2bb = py.numpy.array(xxxxx2b.');

    try
        res = pyrunfile("OptimalExposure.py","q",wwa=wwwa,ww1=www1,ww2=www2,ww0=www0,wwq=wwwq,lamp=lampp,lamn=lamnn,Phip=Phipp,Phin=Phinn,A=AA,xxxxx=xxxxxa,xxxxx2=xxxxx2a,xxxxx2b=xxxxx2bb,alpha=alpha,ky=ky,kz=kz,nlp=nlp,nln=nln,verbose=verbose);
        y=res{1}; yy=double(y);
        z=res{2}; zz=double(z);
        ExposureAcceptability(id).pos=yy;
        ExposureAcceptability(id).mc=[xxxxx zz'];
        fprintf('Day %d solved',id)
    catch ME    
        fprintf('Day %d failed',id)
    end
end

%% Visualization
close all
ndays_job = 10;
for id=1:ndays_job
    try
        yy = ExposureAcceptability(id).pos;
        figure
        plot(xxxxx,A*yy')
        xlim([-0.1,0.1])
        fprintf('Day %d solved\n',id)
        fpath=('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\OptimalExposure\Plots');
        str=strcat('OptimalExposure_OSS2_SPY');
        fname=str;
        saveas(gcf, fullfile(fpath, fname), 'epsc');
    catch ME    
        fprintf('Day %d failed\n',id)
    end
end

%% Deprecated
      % % rebate calculation
      % rhov=(abs(zz).^alpha.*wwa);
      % 
      % % zeta z calculation
      % yz=yyy.*ww1+yyy.*zz.*xxxxx2.*ww2;
      % 
      % obj=sum(yz+rhov)-wwq*yyy;
      % 
      % % constraints
      % conp=ww0'*(max((zz'.*xxxxx2)*ones(1,nlp)-(ones(nxx,1)*(lamp-1)),0))-Phip;
      % conn=ww0'*(max((ones(nxx,1)*(1-lamn)-zz.*xxxxx2)*ones(1,nln)),0))-Phin;
      % 
      % % Bound Constraints 
      % %-100<yy<100
      % % zz>=-1;
       
%% Consider
% DeltaExposure=struct;
% for id=1:ndays
%     disp(id);
%     mat1=DSPEAinputs1(id).probs;
%     avec1=DSPEAinputs1(id).avec;
%     mat2=DSPEAinputs1(id).probs;
%     avec2=DSPEAinputs1(id).avec;
%     xxxxp=mat1(:,1);
%     xxxxn=mat1(:,2);
%     yyap=mat1(:,3);
%     yyan=mat1(:,4);
%     yy1p=mat1(:,5);
%     yy1n=mat1(:,6);
%     yy2p=mat1(:,7);
%     yy2n=mat1(:,8);
%     yy0p=mat1(:,9);
%     yy0n=mat1(:,10);
%     qqp=mat2(:,5);
%     qqn=mat2(:,6);
%     Aap=avec1(1);
%     Aan=avec1(2);
%     A1p=avec1(3);
%     A1n=avec1(4);
%     A2p=avec1(5);
%     A2n=avec1(6);
%     A0p=avec1(7);
%     A0n=avec1(8);
%     AA1p=avec2(5);
%     AA1n=avec2(6);
%     nxp=length(xxxxp);
%     nxn=length(xxxxn);
%     nxx=nxp+nxn;
%     yy=xx;
%     xxxxx=[xxxxn' xxxxp']';
%     A=zeros(nxx,nx);
%     ee=eye(nx);
%     for ic=1:nx
%         A(:,ic)=interp1(xx,ee(:,ic),xxxxx,'spline','extrap');
%     end
%     zz=ones(1,nxx)+.1*randn(1,nxx);
%     zp=max(zz-1,0);
%     zn=max(1-zz,0);
% 
%     yyy=A*yy;
% 
%     wwa=[Aan*yyan' Aap*yyap']';
%     ww1=[A1n*yy1n' A1p*yy1p']';
%     ww2=[A2n*yy2n' A2p*yy2p']';
%     ww0=[A0n*yy0n' A0p*yy0p']';
%     wwq=[AA1n*qqn' AA1p*qqp'];
%     xxxxx=[xxxxn.' xxxxp.'].';
%     xxxxx2=xxxxx.^2;
%     sssss=exp(xxxxx)-1;
%     rebate calculation
%     rhov=((zp'+zn').^alpha.*wwa);
%     
%     % zeta z calculation
%     % yz=yyy.*ww1+yyy.*(zp'-zn').*ww2;
%     %
%     % obj=sum(yz+rhov)-wwq*yyy;
%     
%     % constraints
%     conp=ww0'*(max(((zp'-zn').*xxxxx.^2)*ones(1,nlp)-(ones(nxx,1)*(lamp-1)),0))-Phip;
%     conn=ww0'*(max((ones(nxx,1)*(1-lamn)-((zp'-zn').*xxxxx.^2)*ones(1,nln)),0))-Phin;
%     
%     
%     yz=@(p) p*exp(xxxxx-1).*ww1+p*exp(xxxxx-1).*(zp'-zn').*ww2;
%     obj=-sum(yz(p)+rhov);
%     
%     % optss=optimset('fmincon');
%     % optionss=optimset(optss,'TolX',1e-4,'TolFun',1e-4,'Display','off','MaxFunEvals',3000);
%     ky=1;
%     kzp=1000;
%     kzn=1000;
%     xxxxx2b=1./xxxxx2;
%     pyenv('Version',...
%         'C:\OSSDSP\maximin2\python.exe',...
%         'ExecutionMode','OutOfProcess');
%     py.importlib.import_module('numpy');
% 
%     wwwa = py.numpy.array(wwa.');
%     www1 = py.numpy.array(ww1.');
%     www2 = py.numpy.array(ww2.');
%     www0 = py.numpy.array(ww0.');
%     lampp = py.numpy.array(lamp.');
%     lamnn = py.numpy.array(lamn.');
%     Phipp = py.numpy.array(Phip.');
%     Phinn = py.numpy.array(Phin.');
%     sssssa = py.numpy.array(sssss.');
%     xxxxx2a= py.numpy.array(xxxxx2.');
%     xxxxx2bb = py.numpy.array(xxxxx2b.');
% 
%     verbose='True';
% 
%     try
%         res = pyrunfile("OptimalDeltaExposure.py","z",wwa=wwwa,ww1=www1,ww2=www2,ww0=www0,lamp=lampp,lamn=lamnn,Phip=Phipp,Phin=Phinn,sssss=sssssa,xxxxx2=xxxxx2a,xxxxx2b=xxxxx2bb, alpha=alpha,ky=ky,kzp=kzp,kzn=kzn,nlp=nlp,nln=nln,verbose=verbose);
%         y=res{1}; yy=double(y);
%         zp=res{2}; zzp=double(zp);
%         zn=res{3}; zzn=double(zn);
%         DeltaExposure(id).pos=yy;
%         DeltaExposure(id).mc=[xxxxx zzp' zzn'];
%     catch ME    
%         y=res{1}; yy=double(y);
%         zp=res{2}; zzp=double(zp);
%         zn=res{3}; zzn=double(zn);
%         DeltaExposure(id).pos=yy;
%         DeltaExposure(id).mc=[xxxxx zzp' zzn'];
%     end
% end
% 
% % DeltaExposure Results.
% 
% pvec=zeros(ndays,1);
% for id=1:ndays
%       p=DeltaExposure(id).pos;
%       if isempty(p)==0
%          pvec(id)=p;
%       end
% end
% 
% % x^2 normalized expectations
% mpv=1./bpv;
% gpv=1./bnv;
% vvv=cpv.*gamma(2-ypv).*(1./(mpv-1).^(2-ypv)-1./mpv.^(2-ypv))+cnv.*gamma(2-ynv).*(1./(gpv+1).^(2-ynv)-1./gpv.^(2-ynv));
% 
% qv=quantile(vvv,[.01 .05 .1 .25 .5 .75 .9 .95 .99]);
% 
% % Cone embedding condition
% 
% coneq=@(x) cp.*gamma(2-yp).*(1./(M-x-1).^(2-yp)-1./(M-x).^(2-yp))+cn.*gamma(2-yn).*(1./(G+x+1).^(2-yn)-1./(G+x).^(2-yn));
% etav=zeros(100,1);
% for ic=1:100
%       disp(ic);
%       cp=cpv(ic); M=mpv(ic); yp=ypv(ic);
%       cn=cnv(ic); G=gpv(ic); yn=ynv(ic);
%       coneq=@(x) (cp.*gamma(2-yp).*(1./(M+x-1).^(2-yp)-1./(M+x).^(2-yp))+cn.*gamma(2-yn).*(1./(G+x+1).^(2-yn)-1./(G+x).^(2-yn)))*10000;
%       xx=fsolve(coneq,0,optimset('Display','off'));
%       etav(ic)=xx;
% end
% qeta=quantile(etav,[.01 .05 .1 .25 .5 .75 .9 .95 .99]);