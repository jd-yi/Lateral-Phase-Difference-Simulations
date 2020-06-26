%% 
%This file runs the code to generate Fig 3C. It plots a boxplot that
%calculates the proportions of trials that the attractor network is in
%State 1, where the firing rate of E3>E4 (refer to Materials and Methods
%for details.


%% initialize
clear; 
close all;
%simulation parameters
dt = 0.1;%ms
t = 0:dt:1200;%ms
fs = 1000/dt;
L = length(t);
f = fs*(0:(L/2))/L;
fdiff = diff(f);
fstep = fdiff(1);

%time constants
tau_e = 6;
tau_i = 15;
%refractory constants
r_e = 2/1000;
r_i = 1/1000;

% excitatory part of hypercolum
A_e = zeros(4,length(t));
dA_e = zeros(4,length(t));

% inhibitory part of hypercolum
A_i = zeros(3,length(t));
dA_i = zeros(3,length(t));

%number of reps
nTrials = 100;
%number of phases to test
phases = (pi/180).*[-90 -60 -30 0 30 60 90];
meanof = zeros(nTrials,length(phases));

for kk = 1:length(phases)
for u = 1:nTrials
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
w_in_reccur =-kii.*diag(ones(1,3));

%generate the LFPs
signal = zeros(3,length(t));
fLow = 0.1;
fHigh = 200;
for ff = fLow:fstep:fHigh 
    phase1 = pi*(2*rand()-1);
    shift = phases(kk);%this is the phase difference DELTA THETA (see materials and methods)
    phase2 = phase1 - (shift+0.3.*pi*(2*rand()-1));
    signal(1,:) = signal(1,:) + map(ff,0.5,8,10)*sin(2*pi*(ff/1000)*t-phase1);
    signal(2,:) = signal(2,:) + map(ff,0.5,8,10)*sin(2*pi*(ff/1000)*t-phase2);
end

% noise sources
noiseE = 0.3*randn(4,length(t));
noiseI = 0.3*randn(3,length(t));

% oscillator sources
ampl = 1;%0.3;%0.05;
oscs = zeros(3,length(t));

oscs(1,:)= ampl.*(signal(1,:)./max(abs(signal(1,:))));
oscs(2,:)= ampl.*(signal(2,:)./max(abs(signal(2,:))));

%this code decides the strength of the stimulus. adjust the coefficent for
%different values
stimA = zeros(4,length(t));
stimA(1,:) = 0.5*(t>300&t<1200).*ones(1,length(t));
stimA(2,:) = 0.5*(t>300&t<1200).*ones(1,length(t));

%% simulate
for i = 2:length(t)
    A_e(:,i) = A_e(:,i-1) + dA_e(:,i-1)*dt;
    A_i(:,i) = A_i(:,i-1) + dA_i(:,i-1)*dt;
    
    dA_e(:,i) = (-A_e(:,i-1) + (1-r_e.*A_e(:,i-1)).*max(0,stimA(:,i-1)+w_ex_reccur*A_e(:,i-1)+w_EI_conn*A_i(:,i-1)+noiseE(:,i-1)))/tau_e;
    dA_i(:,i) = (-A_i(:,i-1) + (1-r_i.*A_i(:,i-1)).*max(0,oscs(:,i-1)+w_IE_conn*A_e(:,i-1)+w_in_reccur*A_i(:,i-1)+noiseI(:,i-1)))/tau_i; 
end

%% plot
% 
% figure()
% tiledlayout(4,2)
% nexttile([2,2])
% plot(A_e(3,:),A_e(4,:),'k')
% xlabel('Area 3');
% ylabel('Area 4');
% nexttile([1,2])
% plot(t,A_e(3,:),'r')
% hold on
% plot(t,A_e(4,:),'b')
% legend('Area 3','Area 4','Location','bestoutside')
% xlabel('Time (ms)')
% ylabel('Firing Rate (a.u.)')
% nexttile([1,2])
% plot(t,A_e(3,:)>A_e(4,:),'r')
% xlabel('Time (ms)')
% ylabel('State')
meanof(u,kk)=mean(A_e(3,t>300)>A_e(4,t>300));
end
end
boxplot(meanof,phases.*(180/pi),'PlotStyle','compact','Colors','k')
xlabel('Phase offset (^o)')
ylabel('Proprtion when A_3 > A_4')

%% function defitions 
function P = map(x,sig,c,amp)
    g = amp*(1/2*sqrt(pi*sig))*exp(-0.5*((x-c)/sig).^2);
    P = (1./x.^0.4).*(0.3)+g;
end