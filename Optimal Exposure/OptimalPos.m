clear
clc
%close all

%% Data

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

TT = length(a);

%% Path to OptimalPos_fun file

pythonpath = "OptimalPos_fun_v2.py";

%% Constraints
Cu = 500;
Cl = 500;
lamp = linspace(1,2,Cu);
% lamn1 = linspace(0,1,Cl);
% lamn1 = (1-(1-lamn1).^(1+0.75))/2;
% lamn2 = linspace(0,1,Cl);
% lamn2 = 0.5+(lamn2).^(1+0.75)/2;
% lamn = [lamn1,lamn2];
lamn = linspace(0,1,Cl);
a = 2;
b = 1;
c = 0.5;
gam = 0.75;
Phi_u = zeros(Cu,1);
for i=1:Cu
    Phi_u(i) = Phiup(a,c,gam,lamp(i));
end
Phi_l = Phitil(b,c,lamn);
eps = 0.05;

% Rebate
alpha = 1.2;

%% Discretization 

K = 50; % discretization of y
N = 250; % discretization of z
X = [-2,2];
x = linspace(X(1),X(2),N);
delta = (X(2)-X(1))/N;
x2 = x.*x;
x2inv = 1./x2;

xx = linspace(x(1),x(end),K);
I = eye(K);
MM = zeros(N,K);
for j=1:K
    ej = I(:,j);
    MM(:,j) = interp1(xx,ej,x,'spline','extrap');
end

%% Optimization

y = cell(TT/2);
z = cell(TT/2);

for tt = 1%:TT/2
    % BCGMY
    params = [parmm(2*tt,:),parmm(2*tt+1,:)];
    
    cp = params(1);
    M = 1/params(2);
    yp = params(3);
    cn = params(4);
    G = 1/params(5);
    yn = params(6);
    
    xp = x(x>0);
    xn = x(x<0);
    
%     eta = [0,0];
%     fun = @(eta)NA(cp,M,eta(1),yp,cn,G,eta(2),yn,xp,xn,x,delta);
%     eta = fminunc(fun,eta);
%     M_eta = M+eta(1);
%     G_eta = G+eta(2); %ensuring long stock and inverse stock cannot be scaled up indefinitely.
    
    pa2 = [(cn).*((-xn).^(2*alpha-yn-1)).*exp(G*xn)*delta,...
           (cp).*(xp.^(2*alpha-yp-1)).*exp(-M*xp)*delta];
    p2 = [(cn).*((-xn).^(2-yn-1)).*exp(G*xn)*delta,...
          (cp).*(xp.^(2-yp-1)).*exp(-M*xp)*delta];
    p4 = [(cn).*((-xn).^(4-yn-1)).*exp(G*xn)*delta,...
          (cp).*(xp.^(4-yp-1)).*exp(-M*xp)*delta];
    p0 = [(cn).*((-xn).^(-yn-1)).*exp(G*xn)*delta,...
          (cp).*(xp.^(-yp-1)).*exp(-M*xp)*delta];
    p0inv2 = p0/(p0*x2');
    
    % cost
    cpq = params(7);
    Mq = 1/params(8);
    ypq = params(9);
    cnq = params(10);
    Gq = 1/params(11);
    ynq = params(12);

%     eta = [0,0];
%     fun = @(eta)NA(cpq,Mq,eta(1),ypq,cnq,Gq,eta(2),ynq,xp,xn,x,delta);
%     eta = fminunc(fun,eta);
%     Mqt = Mq+eta(1);
%     Gqt = Gq+eta(2);
%     Mqt = M_eta;
%     Gqt = G_eta;
    
    q0 = [(cnq).*((-xn).^(-ynq-1)).*exp(Gq*xn)*delta,...
          (cpq).*(xp.^(-ypq-1)).*exp(-Mq*xp)*delta]; %Estimated price to unwound the position in two weeks. 
    q0inv2 = q0/(q0*x2');

    % load python.exe
    pyenv('Version',...
        'C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\maximin2\python.exe',...
        'ExecutionMode','OutOfProcess');
    py.importlib.import_module('numpy');
    
    pa2_py = py.numpy.array(pa2.');
    p2_py = py.numpy.array(p2.');
    p4_py = py.numpy.array(p4.');
    p0_py = py.numpy.array(p0.');
    p0inv2_py = py.numpy.array(p0inv2.');
    q0_py = py.numpy.array(q0.');
    MM_py = py.numpy.array(MM.');
    lamp_py = py.numpy.array(lamp.');
    lamn_py = py.numpy.array(lamn.');
    Phi_u_py = py.numpy.array(Phi_u.');
    Phi_l_py = py.numpy.array(Phi_l.');
    x_py = py.numpy.array(x.');
    x2_py = py.numpy.array(x2.');
    x2inv_py = py.numpy.array(x2inv.');

    res = pyrunfile(pythonpath,"z",...
        pa2 = pa2_py,...
        p2 = p2_py,...
        p4 = p4_py,...
        p0 = p0_py,...
        p0inv2 = p0inv2_py,...
        q0 = q0_py,...
        MM = MM_py,...
        lamp = lamp_py,...
        lamn= lamn_py,...
        Phi_u = Phi_u_py,...
        Phi_l = Phi_l_py,...
        x = x_py,...
        x2 = x2_py,...
        x2inv = x2inv_py,...
        N = N,...
        K = K,...
        Cu = Cu,...
        Cl = Cl,...
        alpha = alpha,...
        verbose = 'True');

        y{tt} = MM*transpose(double(res{1}));
        z{tt} = double(res{2});

%     try
%         res = pyrunfile("OptimalPos_fun.py","z",...
%             pa2 = pa2_py,...
%             p2 = p2_py,...
%             p4 = p4_py,...
%             p0 = p0_py,...
%             q0 = q0_py,...
%             MM = MM_py,...
%             lamp = lamp_py,...
%             lamn= lamn_py,...
%             Phi_u = Phi_u_py,...
%             Phi_l = Phi_l_py,...
%             x = x_py,...
%             x2 = x2_py,...
%             x2inv = x2inv_py,...
%             N = N,...
%             K = K,...
%             C = C,...
%             alpha = alpha,...
%             verbose = 'False');
%             res = MM*transpose(double(res));
%     catch
%         fprintf('For tt = %d, the solver CLARABEL failed\n', tt)
%         res = zeros(size(N,1));
%     end

end

%% Visualization

% Optimal position
figure()
% X = [-0.3,0.3];
plot(x,y{tt})
xlim([-0.1 0.1]);
fpath=('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Figures');
str=strcat('OptimalSolution_SPY');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');

figure()
% X = [-0.3,0.3];
plot(x,1+z{tt}.*x2)
xlim([-0.1 0.1]);
fpath=('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Figures');
str=strcat('OptimalSolution_SPY');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');

% Constraints
% zz = z.value
% N = len(a)
% const = []
% for i in range(N): #range(len(a)):
%     const.append(p @ np.maximum(zz-a[i],0))
% 
% const = np.array(const)
% fig = plt.figure()
% axes = fig.add_axes([0.1, 0.1, 0.75, 0.75])
% axes.set_xlim(a[0], a[-1])
% axes.set_ylim(min(const), max(Phi))
% axes.plot(a,const)
% axes.plot(a,Phi)
% # axes.plot(a,dist.Psi(a))
% plt.show()


%% Routines

function phi = Phitil(b,c,lam)
    lamp = lam(lam<1);
    phi(lam<1) = -(1-lamp).*(log((1-lamp)/b)/c)-(b/c)*(1 - (1-lamp)/b );
    phi(lam==1) = -b/c;
    phi(lam==0) = (log(b)+1-b)/c;
end

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

function p = InverseContract(cp,M,yp,cn,G,yn,xp,xn,x,delta)
    nu = [cn*exp(G*xn).*(-xn).^(-1-yn),cp*exp(-M*xp).*xp.^(-1-yp)]*delta;
    p = nu*(exp(-x')-1);
end

function c = NA(cp,M,etap,yp,cn,G,etan,yn,xp,xn,x,delta)
    p = InverseContract(cp,M,yp,cn,G,yn,xp,xn,x,delta);
    M = M+etap;
    G = G+etan;
    nu = [cn*exp(G*xn).*(-xn).^(1-yn),cp*exp(-M*xp).*xp.^(1-yp)]*delta;    
    c = (nu*(exp(-x')-p)).^2+(nu*(exp(x')-1)).^2;
end