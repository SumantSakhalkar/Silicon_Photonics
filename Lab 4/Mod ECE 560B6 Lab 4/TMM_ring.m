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
neff_data=[2.4992 2.4725 2.4449 2.4163, 2.3868];% calculated using MODE
lambda_data=[1500:25:1600]*1e-9;% wavelength points in the MODE simulation
% p = polyfit(lambda_data,neff_data,3);
% neff=p(1)*lambda.^3+p(2)*lambda.^2+p(3)*lambda+p(4);
neff=interp1(lambda_data,neff_data,lambda);
%plot(lambda*1e9, neff)

load('TMM_3DFDTD.mat');
load('TMM_MODE.mat');

%
% Coupler 1
k1=0.223314129; % cross coupling coefficient from Lab 1
Tc1=0.974567596; % coupling loss from Lab 1
t1=sqrt(Tc1^2-k1^2); % straight-through coupling coefficient
C1=-1i/k1*[-t1, 1; -Tc1, t1]; % coupler matrix

%
%Coupler 2
k2=0.22314129; % cross coupling coefficient from Lab 1
Tc2=Tc1;
t2=sqrt(Tc2^2-k2^2);
C2=-1i/k2*[-t2, 1; -Tc2, t2];

%
% Ring cavity
R=9.75e-6; % radius
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

%-------------------------------------------------------------------------------------------------
figure (1); 
box on;
hold on

plot(lambda*1e9, 10*log10(Drop),'r');
legend('MATLAB code')

set(gca, 'Fontsize', 14, 'xlim', [min(lambda)*1e9 max(lambda)*1e9]);
xlabel('wavelength (nm)', 'Fontsize',14); ylabel('Power transmission (dB)', 'Fontsize', 14);

ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';


ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Fontsize', 14,...
    'YAxisLocation','right',...
    'Color','none');
xlabel('wavelength (nm)', 'Fontsize',14); ylabel('Power transmission (dB)', 'Fontsize', 14);

hold on

plot(lambda_MODE*1e9, Tdrop_MODE,'b');
plot(lambda_3DFDTD*1e9, Tdrop_3DFDTD,'g');

legend('MODE','FDTD')
title('Drop Port')

hold off

%--------------------------------------------------------------------------------------------------------

figure (2); 
box on;
hold on

plot(lambda*1e9, 10*log10(Thru),'r');

legend('MATLAB code')
set(gca, 'Fontsize', 14, 'xlim', [min(lambda)*1e9 max(lambda)*1e9]);
xlabel('wavelength (nm)', 'Fontsize',14); ylabel('Power transmission (dB)', 'Fontsize', 14);

ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';


ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'Fontsize', 14,...
    'YAxisLocation','right',...
    'Color','none');
xlabel('wavelength (nm)', 'Fontsize',14); ylabel('Power transmission (dB)', 'Fontsize', 14);

hold on

plot(lambda_MODE*1e9, Tthrough_MODE,'b');
plot(lambda_3DFDTD*1e9, Tthrough_3DFDTD,'g');

legend('MODE','FDTD')


hold off

%--------------------------------------------------------------------------------------------------------
