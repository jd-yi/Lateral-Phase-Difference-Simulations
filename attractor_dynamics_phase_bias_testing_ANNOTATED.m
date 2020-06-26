%% 
%This file runs the code to generate the simulations in Fig 3A,B. It plots
%the phase plane of E3 and E4, time series of E3 and E4, and the State
%time series when E3>E4.

%% initialize
clear; 
close all;
dt = 0.1;
t = 0:dt:2000;
fs = 1000/dt;
L = length(t);
f = fs*(0:(L/2))/L;
fdiff = diff(f);
fstep = fdiff(1);

tau_e = 6;
tau_i = 15;
r_e = 2/1000;
r_i = 1/1000;
% excitatory part of hypercolum
A_e = zeros(4,length(t));
dA_e = zeros(4,length(t));

% inhibitory part of hypercolum
A_i = zeros(3,length(t));
dA_i = zeros(3,length(t));


%% connectivity
% E to I 
kie = 3.5;
w_IE_conn = kie.*[1 0 0 0; 
                  0 1 0 0;
                  0 0 0.5 0.5];

% I to E
kei = 3.25;
w_EI_conn = -kei.*[1 0 0; 
                   0 1 0; 
                   0 0 1; 
                   0 0 1];

% recurrent excitation
kee = 1.5;%1.5;
w_ex_reccur = kee.*(diag(ones(1,4))+[0 0.2 0 0; 
                                     0.2 0 0 0; 
                                     1 0 0 0; 
                                     0 1 0 0]);%0.2

% recurrent inhibition
kii = 2.5;
w_in_reccur =-kii.*[1 0 0; 
                    0 1 0; 
                    0 0 1];

% noise sources
noiseE = 0.3*randn(4,length(t));
noiseI = 0.3*randn(3,length(t));

signal = zeros(3,length(t));

fLow = 0.1;
fHigh = 200;
for ff = fLow:fstep:fHigh 
    phase1 = pi*(2*rand()-1);
    shift = pi/2;%this is the phase difference DELTA THETA (see materials and methods)
    phase2 = phase1 - (shift+0.3.*pi*(2*rand()-1));
    %signal = signal + (1/ff)*sin(2*pi*(ff/1000)*t-phase);
    signal(1,:) = signal(1,:) + map(ff,0.5,8,10)*sin(2*pi*(ff/1000)*t-phase1);
    signal(2,:) = signal(2,:) + map(ff,0.5,8,10)*sin(2*pi*(ff/1000)*t-phase2);
end

% oscillator sources
%lambda = linspace(0,1,400);
ampl = 0.1;%0.05;
oscs = zeros(3,length(t));

oscs(1,:)= ampl.*(signal(1,:)./max(abs(signal(1,:))));
oscs(2,:)= ampl.*(signal(2,:)./max(abs(signal(2,:))));

%this code decides the strength of the stimulus. adjust the coefficent for
%different values
stimA = zeros(4,length(t));
stimA(1,:) = 0.5*(t>500&t<2000).*ones(1,length(t));
stimA(2,:) = 0.5*(t>500&t<2000).*ones(1,length(t));

%% simulate
for i = 2:length(t)
    A_e(:,i) = A_e(:,i-1) + dA_e(:,i-1)*dt;
    A_i(:,i) = A_i(:,i-1) + dA_i(:,i-1)*dt;
    
    dA_e(:,i) = (-A_e(:,i-1) + (1-r_e.*A_e(:,i-1)).*max(0,stimA(:,i-1)+w_ex_reccur*A_e(:,i-1)+w_EI_conn*A_i(:,i-1)+noiseE(:,i-1)))/tau_e;
    dA_i(:,i) = (-A_i(:,i-1) + (1-r_i.*A_i(:,i-1)).*max(0,oscs(:,i-1)+w_IE_conn*A_e(:,i-1)+w_in_reccur*A_i(:,i-1)+noiseI(:,i-1)))/tau_i; 
end
stateVec = A_e(3,:)>A_e(4,:);
%% plot
max3 = max(A_e(3,:));
max4 = max(A_e(3,:));
maxtot = max(max3,max4);
figure()
tiledlayout(4,2)
nexttile([2,2])
pbaspect([1 1 1])
hold on
%phase plane plot
plot(A_e(3,:)./maxtot,A_e(4,:)./maxtot,'k')
plot(0:0.01:1,0:0.01:1,'--r')
hold off
xlabel('Area 3');
ylabel('Area 4');
nexttile([1,2])
%timeseries plot
plot(t,A_e(3,:)./maxtot,'r')
hold on
plot(t,A_e(4,:)./maxtot,'b')
legend('Area 3','Area 4','Location','bestoutside')
xlabel('Time (ms)')
ylabel('Firing Rate (a.u.)')
nexttile([1,2])

%state E3>E4 plot plot
plot(t,stateVec,'r')
xlabel('Time (ms)')
ylabel('State')

%% 3d plotting (not in manuscript!)
figure()
plot3(A_e(3,:),A_e(4,:),A_i(3,:),'r','LineWidth',1)
xlabel('Area 3')
ylabel('Area 4')
zlabel('Inhibitory Pool')

%% function defitions 
function P = map(x,sig,c,amp)
    g = amp*(1/2*sqrt(pi*sig))*exp(-0.5*((x-c)/sig).^2);
    P = (1./x.^0.4).*(0.3)+g;
    %P = g;
end