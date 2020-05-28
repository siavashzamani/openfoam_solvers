%Omega Analysis
clear all; clc;

for i=1:1:10
ID = num2str(i);
omega_ID = '1';
filename = strcat('omega_',omega_ID,'.',ID,'.csv');
A = csvread(filename,1);
U_max(i) = max(A(:,1));
time(i) = (i)*2*0.001;
end

U_max_mean = mean(U_max);