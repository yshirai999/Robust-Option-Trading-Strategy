clear
clc
close all

load('C:\Users\yoshi\OneDrive\Desktop\Research\OptimalDerivativePos\Maximin\Data\OSS2RNSOMoutput.mat')

n = length(ODPBG(1).dens1(:,1));
p = cell(50);
q = cell(50);
for i = 1:50
    p{i} = ODPBG(i).dens1;
    q{i} = ODPBG(i).dens2;
end

test = 1;
figure
hold on
plot(p{test}(:,1),p{test}(:,2))
plot(q{test}(:,1),q{test}(:,2))
legend('$p$','$q$','interpreter','latex')