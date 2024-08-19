clear
clc

%% Parameters

% BG parameters
params = [0.04,10/52,1,5/52,0.02,10/52,1.5,5/52];
bp = params(1);
cp = params(2);
bn = params(3);
cn = params(4);

% Discretization
k = 8;
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
k = 8;
M = [-0.5,0.5];

b = (1/bptil+1/bntil);
gammap = gamma(cptil);
gamman = gamma(cntil);
lam = 0.5*(cptil+cntil);
mu = 0.5*(cptil+cntil-1);
x = y*b;
xp = x(ip);
xm = -x(im);

W = zeros(size(x));
p = zeros(size(x));
W(ip) = exp(-0.5*xp) .* (xp.^(mu+0.5)) .* kummerU(mu-lam+0.5, 1+2*mu, xp);
W(im) = exp(-0.5*xm) .* (xm.^(mu+0.5)) .* kummerU(mu-lam+0.5, 1+2*mu, xm);
q(ip) = ( (bptil)^(-cptil) ) * ( (bntil)^(-cntil) ) * ( (yp).^(0.5*(cptil+cntil)-1) ) .* exp(-0.5*yp) .* W(ip) / gammap;
q(im) = ( (bptil)^(-cptil) ) * ( (bntil)^(-cntil) ) * ( (ym).^(0.5*(cptil+cntil)-1) ) .* exp(-0.5*ym) .* W(im) / gamman;
q = q * delta / ( b^(0.5*(cptil+cntil)) );
q = q/sum(p);

% Constraints
N = 50;
a = linspace(0.5,2,N);
lam = 0.25;
Phi = MMV_Phi(a,lam);

% Rebate
theta = 0.25;
P = diag(p);

%% Optimization

% load python libraries

pyenv('Version',...
    'C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\maximin2\python.exe',...
    'ExecutionMode','OutOfProcess');
py.importlib.import_module('dsp');
py.importlib.import_module('cvxpy');
py.importlib.import_module('numpy');

l_mat = rand(2,2);
l_py = py.numpy.random.random([int64(2), int64(2)]);

npA = py.numpy.array(l_mat(:).');

p = py.numpy.array(p.');
q = py.numpy.array(q.');
x = py.numpy.array(y.');

P = py.numpy.diag(p);

pyrun("f = dsp.inner(z, P @ (y - q @ y))", " f_mat");
% rho = p @ (cp.multiply(theta, cp.power(z,alpha))+cp.multiply(1-theta,cp.power(z,-beta)))
% obj = dsp.MinimizeMaximize(rho+f)
% constraints = [p @ z == 1, z >= 0]
% for i in range(N): #range(len(a)):
%     constraints.append(p @ cp.maximum(z-a[i],0) <= Phi[i])

%pyrun("l = py.numpy.array(l)")

%% Routines

function psi = MMV_Psi(x,lam)
    psi = 1 - ( 1-x.^(1/(lam+1)) ).^(1+lam);
end

function phi = MMV_Phi(a,lam)
    phi = ( 1 ./ a.^(1/lam) ).^(lam+1);
end