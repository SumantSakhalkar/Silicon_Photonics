% Microring resonator model using the transfer matrix method (TMM)
%
% by Wei Shi -- wei.j.shi@gmail.com
% created in 2012 and modified in May 2013
%
clear 
close all
%
lambda=(1500:0.01:1600)*1e-9;
%lambda=(1550:0.001:1552)*1e-9;
%
% neff_data=[2.408346 2.377321 2.34625 2.315149, 2.284034];%(450,220) calculated using MODE 
% neff_data=[2.427555 2.397263 2.366924 2.33655, 2.306156];%(460,220) calculated using MODE 
% neff_data=[2.390731 2.359011 2.327248 2.295459, 2.263668];%(440,220) calculated using MODE
% neff_data=[2.426911 2.396094 2.365215 2.334288, 2.303329];%(450,225) calculated using MODE
% neff_data=[2.446166 2.41609 2.385954 2.355766, 2.32554];%(460,225) calculated using MODE
% neff_data=[2.409996 2.378502 2.346946 2.315345, 2.283718];%(440,225) calculated using MODE
% neff_data=[2.390037 2.358805 2.327507 2.296196, 2.264894];%(450,215) calculated using MODE
% neff_data=[2.409487 2.378953 2.348387 2.317804, 2.287218];%(460,215) calculated using MODE
neff_data=[2.372962 2.34102 2.309052 2.277079, 2.245125];%(440,215) calculated using MODE

lambda_data=[1500:25:1600]*1e-9;% wavelength points in the MODE simulation
% p = polyfit(lambda_data,neff_data,3);
% neff=p(1)*lambda.^3+p(2)*lambda.^2+p(3)*lambda+p(4);
neff=interp1(lambda_data,neff_data,lambda);
%plot(lambda*1e9, neff)

% load('TMM_3DFDTD.mat');
% 
%
% Coupler 1
% k1=0.0101918; % (450,220)cross coupling coefficient from Lab 1
% k1=0.00959364; % (460,220)cross coupling coefficient from Lab 1
% k1=0.0107178; % (440,220)cross coupling coefficient from Lab 1
% k1=0.00951921; % (450,225)cross coupling coefficient from Lab 1
% k1=0.00897301; % (460,225)cross coupling coefficient from Lab 1
% k1=0.00999707; % (440,225)cross coupling coefficient from Lab 1
% k1=0.0108501; % (450,215)cross coupling coefficient from Lab 1
% k1=0.0101981; % (460,215)cross coupling coefficient from Lab 1
k1=0.0114195; % (440,215)cross coupling coefficient from Lab 1

% Tc1=0.989506; % (450,220)coupling loss from Lab 1
% Tc1=0.990189; % (460,220)coupling loss from Lab 1
% Tc1=0.988954; % (440,220)coupling loss from Lab 1
% Tc1=0.99022; % (450,225)coupling loss from Lab 1
% Tc1=0.990848; % (460,225)coupling loss from Lab 1
% Tc1=0.989679; % (440,225)coupling loss from Lab 1
% Tc1=0.988822; % (450,215)coupling loss from Lab 1
% Tc1=0.98954; % (460,215)coupling loss from Lab 1
Tc1=0.988247; % (440,215)coupling loss from Lab 1

t1=sqrt(Tc1^2-k1^2); % straight-through coupling coefficient
C1=-1i/k1*[-t1, 1; -Tc1, t1]; % coupler matrix

%
%Coupler 2
k2=k1; % cross coupling coefficient from Lab 1
Tc2=Tc1;
t2=sqrt(Tc2^2-k2^2);
C2=-1i/k2*[-t2, 1; -Tc2, t2];

%
% Ring cavity
R=10e-6; % radius
L_rt=2*pi*R; % roundtrip length
alpha=4; % loss in dB/cm
A=10^(-alpha*L_rt*100/10); % roundtrip loss
a=sqrt(A);
beta=neff*2*pi./lambda; % propogation constant

%
% % Coupler 2 redefined for critical coupling
% t2=t1/a;
% Tc2=Tc1;
% k2=sqrt(Tc2-t2^2);
% C2=-1i/k2*[-t2, 1; -Tc2, t2];
% %

Thru=zeros(1, length(lambda)); % through-port response
Drop=zeros(1, length(lambda)); % drop-port response
%

for ii=1:length(lambda)
    P=[0, sqrt(a)*exp(-1i*beta(ii)*L_rt/2);...
        1/(sqrt(a)*exp(-1i*beta(ii)*L_rt/2)), 0]; % ring cavity matrix
    H=C1*P*C2; % total transfer matrix
    Drop(ii)=abs(1/H(1,2))^2; % drop-port response
    Thru(ii)=abs(H(2,2)/H(1,2))^2; % through-port response
end
save('matlab440215.mat');
