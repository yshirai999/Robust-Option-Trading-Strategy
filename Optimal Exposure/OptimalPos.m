clear
clc
close all

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
%% Constraints
C = 50;
lamp = linspace(1.1,2,C);
lamn = linspace(0.1,0.9,C);
a = 1;
b = 1;
c = 0.01;
gam = 1;
Phi_u = zeros(C,1);
for i=1:C
    Phi_u(i) = Phiup(a,c,gam,lamp(i));
end
Phi_l = Phitil(b,c,lamn);
eps = 0.01;

% Rebate
alpha = 1.2;

%% Discretization 

K = 50; % discretization of y
N = 1000; %discretization of z
X = [-1,1];
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

for tt = 2%1:TT/2
    % BCGMY
    params = [parmm(2*tt-1,:),parmm(2*tt,:)];
    %[0.04,13/52,1.2,2/52,0,0,0.02,10/52,1.5,5/52,0,0];
    
    cp = params(1);
    M = 1/params(2);
    yp = params(3);
    cn = params(4);
    G = 1/params(5);
    yn = params(6);
    
    xp = x(x>0);
    xn = x(x<0);

    pa2 = [(cn).*((-xn).^(2*alpha-yn-1)).*exp(G*xn)*delta,...
           (cp).*(xp.^(2*alpha-yp-1)).*exp(-M*xp)*delta];
    p2 = [(cn).*((-xn).^(2-yn-1)).*exp(G*xn)*delta,...
          (cp).*(xp.^(2-yp-1)).*exp(-M*xp)*delta];
    p4 = [(cn).*((-xn).^(4-yn-1)).*exp(G*xn)*delta,...
          (cp).*(xp.^(4-yp-1)).*exp(-M*xp)*delta];
    p0 = [(cn).*((-xn).^(-yn-1)).*exp(G*xn).*(xn<-eps)*delta,...
          (cp).*(xp.^(-yp-1)).*exp(-M*xp).*(xp>eps)*delta];
    
    % cost
    cpq = params(7);
    Mq = 1/params(8);
    ypq = params(9);
    cnq = params(10);
    Gq = 1/params(11);
    ynq = params(12);
    
    q0 = [(cnq).*((-xn).^(-ynq-1)).*exp(Gq*xn).*(xn<-eps)*delta,...
          (cpq).*(xp.^(-ypq-1)).*exp(-Mq*xp).*(xp>eps)*delta];

    % load python.exe
    pyenv('Version',...
        'C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\maximin2\python.exe',...
        'ExecutionMode','OutOfProcess');
    py.importlib.import_module('numpy');
    
    pa2_py = py.numpy.array(pa2.');
    p2_py = py.numpy.array(p2.');
    p4_py = py.numpy.array(p4.');
    p0_py = py.numpy.array(p0.');
    q0_py = py.numpy.array(q0.');
    MM_py = py.numpy.array(MM.');
    lamp_py = py.numpy.array(lamp.');
    lamn_py = py.numpy.array(lamn.');
    Phi_u_py = py.numpy.array(Phi_u.');
    Phi_l_py = py.numpy.array(Phi_l.');
    x2_py = py.numpy.array(x2.');
    x2inv_py = py.numpy.array(x2inv.');
    
    res = pyrunfile("OptimalPos_fun.py","z",...
        pa2 = pa2_py,...
        p2 = p2_py,...
        p4 = p4_py,...
        p0 = p0_py,...
        q0 = q0_py,...
        MM = MM_py,...
        lamp = lamp_py,...
        lamn= lamn_py,...
        Phi_u = Phi_u_py,...
        Phi_l = Phi_l_py,...
        x2 = x2_py,...
        x2inv = x2inv_py,...
        N = N,...
        K = K,...
        C = C,...
        alpha = alpha,...
        verbose = 'True');
        res = A*transpose(double(res));

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
%             x2 = x2_py,...
%             x2inv = x2inv_py,...
%             N = N,...
%             K = K,...
%             C = C,...
%             alpha = alpha,...
%             verbose = 'False');
%             res = A*transpose(double(res));
%     catch
%         fprintf('For tt = %d, the solver CLARABEL failed\n', tt)
%         res = zeros(size(N,1));
%     end

end

%% Visualization

% % Optimal position
% figure()
% X = [-0.3,0.3];
% xlim = ([log(W)+X(1), log(W)+X(2)]);
% plot(x,res)
% fpath=('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Figures');
% str=strcat('OptimalSolution_SPY');
% fname=str;
% saveas(gcf, fullfile(fpath, fname), 'epsc');

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
    phi = -(1-lam).*(log((1-lam)/b)/c)-(b/c)*(1 - (1-lam)/b );
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