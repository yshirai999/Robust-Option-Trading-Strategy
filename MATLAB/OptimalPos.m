clear
clc
close all

%% Parameters

% Constraints
N = 50;
a = linspace(0.5,2,N);
lam = 0.25;
Phi = MMV_Phi(a,lam);

% Rebate
theta = 0.75;
alpha = 1.25;
beta = 0.25;

%% Densities

load('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\Data\OSS2RNSOMoutput.mat')

n = length(ODPBG(1).dens1(:,1));
T = length(ODPBG);
p = cell(50);
q = cell(50);
for i = 1:50
    p{i} = ODPBG(i).dens1;
    q{i} = ODPBG(i).dens2;
end

%% Optimization

% load python.exe

pyenv('Version',...
    'C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\maximin2\python.exe',...
    'ExecutionMode','OutOfProcess');
py.importlib.import_module('numpy');

%%

y = cell(T);
z = cell(T);
for i = 1:1
    xx = p{i}(:,1);
    pp = p{i}(:,2)/0.0005;
    qq = q{i}(:,2)/0.0005;
    
    % Interpolation matrix
    k = length(xx);
    ky = 40;
    [~,ind] = maxk(pp,ky);
    I = eye(ky);
    xxx = sort(xx(ind));
    A = zeros(k,ky);
    for j=1:ky
        ej = I(:,j);
        A(:,j) = interp1(xxx,ej,xx,'spline','extrap');
    end

    ppp = py.numpy.array(pp.');
    qqq = py.numpy.array(qq.');
    aa = py.numpy.array(a.');
    AA = py.numpy.array(A.');
    Phipy = py.numpy.array(Phi.');
    
    res = pyrunfile("OptimalPos_fun.py","z",p=ppp,q=qqq,k=k,ky=ky,a=aa,A=AA,Phi=Phipy,theta=theta,alpha=alpha,beta=beta,N=N);
    
    yy = double(res{1});
    y{i} = A*yy;
    z{i} = double(res{2});

end



%% Visualization

% Optimal position
figure()
plot(p{i}(:,1),y{1})
fpath=('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Figures');
str=strcat('ODP_SPY');
fname=str;
saveas(gcf, fullfile(fpath, fname), 'epsc');


fprintf('Cost of implementing strategy y is %d\n', q * y{i}')

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