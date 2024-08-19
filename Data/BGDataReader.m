clc
clear

a = [512,1024];
for i = 1:2
    filename = strcat('bgset',num2str(a(i)));
    A = importdata(strcat(filename,'.dat'), ' ', 1);
    save(filename,"A")
end