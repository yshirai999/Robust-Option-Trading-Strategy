clear
clc
close all

% Parameters

% BG parameters
params = [0.04,10/52,1,5/52,0.02,10/52,1.5,5/52];
bp = params(1);
cp = params(2);
bn = params(3);
cn = params(4);

% Discretization
k = 10;
M = [-0.5,0.5];
y = linspace(M(1),M(2),2^k);
ip = (y > 0);
im = (y < 0);
yp = y(ip);
ym = -y(im);
delta = (M(2)-M(1))/(2^k-1);

b = (1/bp+1/bn);
gammap = gamma(cp);
gamman = gamma(cn);
lam = 0.5*(cp+cn);
mu = 0.5*(cp+cn-1);
x = y*b;
xp = x(ip);
xm = -x(im);

W = zeros(size(x));
p = zeros(size(x));
W(ip) = exp(-0.5*xp) .* (xp.^(mu+0.5)) .* kummerU(mu-lam+0.5, 1+2*mu, xp);
W(im) = exp(-0.5*xm) .* (xm.^(mu+0.5)) .* kummerU(mu-lam+0.5, 1+2*mu, xm);
p(ip) = ( (bp)^(-cp) ) * ( (bn)^(-cn) ) * ( (yp).^(0.5*(cp+cn)-1) ) .* exp(-0.5*yp) .* W(ip) / gammap;
p(im) = ( (bp)^(-cp) ) * ( (bn)^(-cn) ) * ( (ym).^(0.5*(cp+cn)-1) ) .* exp(-0.5*ym) .* W(im) / gamman;
p = p * delta / ( b^(0.5*(cp+cn)) );
p = p/sum(p);

% BG Risk neutral parameters
bptil = params(5);
cptil = params(6);
bntil = params(7);
cntil = params(8);

% Discretization
b = (1/bptil+1/bntil);
gammap = gamma(cptil);
gamman = gamma(cntil);
lam = 0.5*(cptil+cntil);
mu = 0.5*(cptil+cntil-1);
x = y*b;
xp = x(ip);
xm = -x(im);

W = zeros(size(x));
q = zeros(size(x));
W(ip) = exp(-0.5*xp) .* (xp.^(mu+0.5)) .* kummerU(mu-lam+0.5, 1+2*mu, xp);
W(im) = exp(-0.5*xm) .* (xm.^(mu+0.5)) .* kummerU(mu-lam+0.5, 1+2*mu, xm);
q(ip) = ( (bptil)^(-cptil) ) * ( (bntil)^(-cntil) ) * ( (yp).^(0.5*(cptil+cntil)-1) ) .* exp(-0.5*yp) .* W(ip) / gammap;
q(im) = ( (bptil)^(-cptil) ) * ( (bntil)^(-cntil) ) * ( (ym).^(0.5*(cptil+cntil)-1) ) .* exp(-0.5*ym) .* W(im) / gamman;
q = q * delta / ( b^(0.5*(cptil+cntil)) );
q = q/sum(q);

% Interpolation matrix
ky = 5;
[~,ind] = maxk(p,2^ky);
I = eye(2^ky);
yy = sort(y(ind));
A = zeros(2^k,2^ky);
for j=1:2^ky
    ej = I(:,j);
    A(:,j) = interp1(yy,ej,y,'spline','extrap');
end


%%
% Constraints
N = 50;
a = linspace(0.5,2,N);
lam = 0.25;
Phi = MMV_Phi(a,lam);

% Rebate
theta = 0.75;
alpha = 1.25;
beta = 0.25;

%% Optimization

% load python.exe

pyenv('Version',...
    'C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\maximin2\python.exe',...
    'ExecutionMode','OutOfProcess');
py.importlib.import_module('numpy');

p = py.numpy.array(p.');
q = py.numpy.array(q.');
x = py.numpy.array(y.');
a = py.numpy.array(a.');
Anp = py.numpy.array(A.');
Phi = py.numpy.array(Phi.');

res = pyrunfile("OptimalPos_fun.py","z",p=p,q=q,x=x,a=a,Phi=Phi,A=Anp,k=k,ky=ky,theta=theta,alpha=alpha,beta=beta,W=W,N=N);

res = A*transpose(double(res));
q = double(q);

%% Visualization

% Optimal position
figure()
M = [-0.3,0.3];
xlim = ([log(W)+M(1), log(W)+M(2)]);
plot(x,res)
fpath=('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Figures');
str=strcat('OptimalSolution_SPY');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');


fprintf('Cost of implementing strategy y is %d\n', q * res)

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

function psi = MMV_Psi(x,lam)
    psi = 1 - ( 1-x.^(1/(lam+1)) ).^(1+lam);
end

function phi = MMV_Phi(a,lam)
    phi = ( 1 ./ a.^(1/lam) ).^(lam+1);
end