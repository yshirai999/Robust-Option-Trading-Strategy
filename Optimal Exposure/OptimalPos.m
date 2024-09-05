clear
clc
close all

%% Constraints
C = 50;
lamp = linspace(1.1,2,C);
lamn = linspace(0.1,0.9,C);
a = 0.01;
b = 1;
c = 0.5;
gam = 1;
Phi_u = zeros(C,1);
for i=1:C
    Phi_u(i) = Phiup(a,c,gam,lamp(i));
end
Phi_l = Phitil(b,c,lamn);
eps = 0.01;

% Rebate
theta = 0.25;
alpha = 1.2;
beta = 0.25;

%% BCGMY
params = [0.04,13/52,1.2,2/52,0,0,0.02,10/52,1.5,5/52,0,0];
M = 1/params(1);
G = 1/params(3);
cp = params(2);
cn = params(4);
yp = params(5);
yn = params(6);

% Discretization
K = 50; % discretization of y
N = 10000; %discretization of z
X = [-0.5,0.5];
x = linspace(X(1),X(2),N);
x2 = x.*x;
x2inv = 1./x2;

% alpha2-rebate
A_pa2 = cp*M^(-2*alpha+yp).*gamma(2*alpha-yp);
p_pa2 = (cp/A_pa2).*(x.^(2*alpha-yp-1)).*exp(-M*x).*(x>0);
A_na2 = cn*G.^(-2*alpha-yn).*gamma(2*alpha-yn);
p_na2 = (cn/A_na2).*((-x).^(2*alpha-yn-1)).*exp(G*x).*(x<0);

% inner2
A_p2 = cp*M.^(-2-yp).*gamma(2-yp);
p_p2 = (cp/A_p2).*(x.^(2-yp-1)).*exp(-M*x).*(x>0);
A_n2 = cn*G.^(-2-yn).*gamma(2-yn);
p_n2 = (cn/A_n2).*((-x).^(2-yn-1)).*exp(G*x).*(x<0);

% inner4
A_p4 = cp*M.^(-4-yp).*gamma(4-yp);
p_p4 = (cp/A_p4).*(x.^(4-yp-1)).*exp(-M*x).*(x>0);
A_n4 = cn*G.^(-4-yn).*gamma(4-yn);
p_n4 = (cn/A_n4).*((-x).^(4-yn-1)).*exp(G*x).*(x<0);

% constraints
if yp>0
    B_p = (cp*M^(-yp)) * ( (exp(-eps)*(eps^(-yp))/yp) - (1/yp)*igamma(yp,eps) );
else
    B_p = (cp*M^(-yp)) * expint(eps);
end
p_p0 = (cp/B_p)*(exp(-M*x)./(x.^(1+yp))).*(x>eps);
if yn>0
    B_n = (cp*G^(-yn)) * ( (exp(-eps)*(eps^(-yn))/yn) - (1/yn)*igamma(yn,eps) );
else
    B_n = (cn*G^(-yn)) * expint(eps);
end
p_n0 = (cn/B_n)*(exp(G*x)./((-x).^(1+yn))).*(x<-eps);

% cost
Mq = 1/params(7);
Gq = 1/params(8);
cpq = params(9);
cnq = params(10);
ypq = params(11);
ynq = params(12);

q_p0 = (cpq)*(exp(-Mq*x)./(x.^(1+ypq))).*(x>eps);
q_n0 = (cnq)*(exp(Gq*x)./((-x).^(1+ynq))).*(x<-eps);

%% Interpolation matrix
xx = linspace(x(1),x(end),K);
I = eye(K);
M = zeros(N,K);
for j=1:K
    ej = I(:,j);
    M(:,j) = interp1(xx,ej,x,'spline','extrap');
end


%% Optimization

% load python.exe

pyenv('Version',...
    'C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\maximin2\python.exe',...
    'ExecutionMode','OutOfProcess');
py.importlib.import_module('numpy');

p_pa2_py = py.numpy.array(p_pa2.');
p_na2_py = py.numpy.array(p_na2.');
p_p2_py = py.numpy.array(p_p2.');
p_n2_py = py.numpy.array(p_n2.');
p_p4_py = py.numpy.array(p_p4.');
p_n4_py = py.numpy.array(p_n4.');
A_pa2_py = py.numpy.array(A_pa2.');
A_na2_py = py.numpy.array(A_na2.');
A_p2_py = py.numpy.array(A_p2.');
A_n2_py = py.numpy.array(A_n2.');
A_p4_py = py.numpy.array(A_p4.');
A_n4_py = py.numpy.array(A_n4.');
B_p_py = py.numpy.array(B_p.');
p_p0_py = py.numpy.array(p_p0.');
B_n_py = py.numpy.array(B_n.');
p_n0_py = py.numpy.array(p_n0.');
q_p0_py = py.numpy.array(q_p0.');
q_n0_py = py.numpy.array(q_n0.');
M_py = py.numpy.array(M.');
lamp_py = py.numpy.array(lamp.');
lamn_py = py.numpy.array(lamn.');
Phi_u_py = py.numpy.array(Phi_u.');
Phi_l_py = py.numpy.array(Phi_l.');
x2_py = py.numpy.array(x2.');
x2inv_py = py.numpy.array(x2inv.');

res = pyrunfile("OptimalPos_fun.py","z",...
    p_pa2 = p_pa2_py,...
    p_na2 = p_na2_py,...
    p_p2 = p_p2_py,...
    p_n2 = p_n2_py,...
    p_p4 = p_p2_py,...
    p_n4 = p_n4_py,...
    A_pa2 = A_pa2_py,...
    A_na2 = A_na2_py,...
    A_p2 = A_p2_py,...
    A_n2 = A_n2_py,...
    A_p4 = A_p4_py,...
    A_n4 = A_n4_py,...
    B_p = B_p_py,...
    p_p0 = p_p0_py,...
    B_n = B_n_py,...
    p_n0 = p_n0_py,...
    q_p0 = q_p0_py,...
    q_n0 = q_n0_py,...
    M = M_py,...
    lamp = lamp_py,...
    lamn= lamn_py,...
    Phi_u = Phi_u_py,...
    Phi_l = Phi_l_py,...
    x2 = x2_py,...
    x2inv = x2inv_py,...
    N = N,...
    K = K,...
    C = C,...
    theta = theta,...
    alpha = alpha);

res = A*transpose(double(res));

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