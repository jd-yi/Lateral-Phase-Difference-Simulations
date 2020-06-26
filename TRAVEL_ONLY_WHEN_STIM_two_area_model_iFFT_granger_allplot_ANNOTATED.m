%% 
%This file runs the code to generate the simulations in Fig. 2. It plots
%the coherency, Granger Causality, DAI, and PAC via MI results (see the
%materials and methods.

% REQUIRES the Wavelet Toolbox, MVGC Toolbox, CHRONUX Toolbox, and the shadedErrorBar function. Please refer
% to manuscript for references and more details. 
%% initialize
clear
dt = 0.1;
t = 0:dt:2000;
fs = 1000/dt;
L = length(t);
f = fs*(0:(L/2))/L;
fdiff = diff(f);
fstep = fdiff(1);

%number of populations to consider (only use 2 here)
N = 2;

tau_e = 6;
tau_i = 15;
r_e = 2/1000;
r_i = 1/1000;
% excitatory part of hypercolum
A_e = zeros(N,length(t));
dA_e = zeros(N,length(t));

% inhibitory part of hypercolum
A_i = zeros(N,length(t));
dA_i = zeros(N,length(t));

% E to I 
kie = 3.5;
w_IE_conn = kie.*(diag(ones(1,N)));

% I to E
kei = 3.25;
w_EI_conn = kei.*(diag(-1.*ones(1,N)));

% recurrent excitation
kee = 1.5;
w_ex_reccur = kee.*(diag(ones(1,N))+0.2.*(eye(N,N)==0));%0.2

% recurrent inhibition
kii = 2.5;
w_in_reccur =kii.*diag(-1.*ones(1,N));

%stores the timeseries for each trial for analysis later
nTrials = 30;
data1Matrix = NaN(length(t(t>100)),nTrials);
data2Matrix = NaN(length(t(t>100)),nTrials);

for f = 1:nTrials
% noise sources
noiseE = 0.3*randn(N,length(t));%400*randn(N,length(t));
noiseI = 0.3*randn(N,length(t));%200*randn(N,length(t));

signal1 = zeros(1,length(t));
signal2 = zeros(1,length(t));

fLow = 0.1;
fHigh = 200;
for ff = fLow:fstep:fHigh 
    phase1 = pi*(2*rand()-1);
    shift = -pi/2;%this parameter determines the phase shift of the two regions
    phase2 = phase1 - (shift+0.3.*pi*(2*rand()-1));
    %signal = signal + (1/ff)*sin(2*pi*(ff/1000)*t-phase);
    signal1 = signal1 + map(ff,0.5,8,10)*sin(2*pi*(ff/1000)*t-phase1);
    signal2 = signal2 + map(ff,0.5,8,10)*sin(2*pi*(ff/1000)*t-phase2);
end

% oscillator sources

%this is the amplitude of the LFP modulations 
ampl = 0.1;%2;%0.05;
oscs = NaN(N,length(t));

oscs(1,:)= ampl.*(signal1./max(abs(signal1)));
oscs(2,:)= ampl.*(signal2./max(abs(signal2)));

%stimulations
tON = 0;
tOFF = 5000;

%adjust the coefficent to change the amplitude of the stimulation provided
%to the two regions. 
isON = (t>=tON&t<=tOFF);
stimA(1,:) = 0.5*isON.*ones(1,length(t));
stimA(2,:) = 0.5*isON.*ones(1,length(t));
%% simulate

for i = 2:length(t)
    A_e(:,i) = A_e(:,i-1) + dA_e(:,i-1)*dt;
    A_i(:,i) = A_i(:,i-1) + dA_i(:,i-1)*dt;
    
    dA_e(:,i) = (-A_e(:,i-1) + (1-r_e.*A_e(:,i-1)).*max(0,stimA(:,i-1)+w_ex_reccur*A_e(:,i-1)+w_EI_conn*A_i(:,i-1)+noiseE(:,i-1)))/tau_e;
    dA_i(:,i) = (-A_i(:,i-1) + (1-r_i.*A_i(:,i-1)).*max(0,oscs(:,i-1)+w_IE_conn*A_e(:,i-1)+w_in_reccur*A_i(:,i-1)+noiseI(:,i-1)))/tau_i;
end
    data1Matrix(:,f)=A_e(1,t>100);%+A_i(1,t>100);
    data2Matrix(:,f)=A_e(2,t>100);%+A_i(2,t>100);
end
data1 = data1Matrix-mean(data1Matrix,1);%detrend(data1Matrix,1);
data2 = data2Matrix-mean(data2Matrix,1);%detrend(data2Matrix,1);

%% coherency (CHRONUX)
params.Fs = fs;
params.err = [2 0.01];
params.pad = 1;
params.trialave = 1;
params.fpass = [0 100];
params.tapers = [3 5];
[C,phi,S12,S1,S2,ff,confC,phistd,Cerr]=coherencyc(data1,data2,params);
shadeUP = Cerr(1,:)-C';
shadeDN = C'-Cerr(2,:);
shadedErrorBar(ff,C,[shadeUP; shadeDN],'lineProps','-k')
hold on 
line([0 100],[confC confC],'Color','k','LineStyle','--','LineWidth',1)
hold off
ylabel('Cohernce')
xlabel('Frequency (Hz)')

%% granger causality (MVGC)
if 1==1%allows you to turn it on and off (it can take a long time to run)
X(1,:,:) = data1;
X(2,:,:) = data2;
momax = 50;
% % RUN THE BIC if you want to know the model order (how many lags to use)
% [~,BIC] = tsdata_to_infocrit(X,30,'LWR');
% [~,bmo_BIC] = min(BIC);
% 
% number of freq bands to consider
freqres = 10000;
freqs = sfreqs(freqres,fs);

%number of bootstrap samples
nSamp = 10;

fBpw = bootstrap_tsdata_to_spwcgc(X,17,freqres,nSamp,0,1E-6);%27

%GRANGER 2 -> 1
fB12pw = reshape(fBpw(:,1,2,:),nSamp,[]);
meanfB12pw = mean(fB12pw);
SEfB12pw = std(fB12pw);
%GRANGER 1 -> 2
fB21pw = reshape(fBpw(:,2,1,:),nSamp,[]);
meanfB21pw = mean(fB21pw);
SEfB21pw = std(fB21pw);

%% plot granger

figure()
tiledlayout(3,2)
nexttile([2 2])
shadedErrorBar(freqs(freqs<=100),meanfB12pw(freqs<=100),SEfB12pw(freqs<=100),'lineprops','-b');
hold on
shadedErrorBar(freqs(freqs<=100),meanfB21pw(freqs<=100),SEfB21pw(freqs<=100),'lineprops','-r');
xlabel('Frequency (Hz)');
ylabel('Granger Causality')
legend('2->1','1->2')

Hfreq = ttest(fB12pw(:,freqs<100),fB21pw(:,freqs<100),0.01);

%Directed Asymmetry Index  1 to 2
DAI = mean((fB21pw-fB12pw)./(fB21pw+fB12pw));
DAISE = std((fB21pw-fB12pw)./(fB21pw+fB12pw));
nexttile([1 2])
shadedErrorBar(freqs(freqs<=100),DAI(freqs<=100),DAISE(freqs<=100),'lineprops','-k');
xlabel('Frequency (Hz)');
ylabel('DAI 1->2')
line([0 100],[0 0],'Color','k','LineStyle','--')
end

%% cwt analysis
tintest = 100;
fintrest = 10;
figure;
%[wtsg1,fsg1]=cwt(mean(data1,2),fs,'FrequencyLimits', [0 100],'VoicesPerOctave',48);
%[wtsg2,fsg2]=cwt(mean(data2,2),fs,'FrequencyLimits', [0 100],'VoicesPerOctave',48);
[wtsg1,fsg1]=cwt(data1(:,30),fs,'FrequencyLimits', [0 100],'VoicesPerOctave',48);
[wtsg2,fsg2]=cwt(data2(:,30),fs,'FrequencyLimits', [0 100],'VoicesPerOctave',48);
tcwt = t((t>tintest));

tiledlayout(2,1)
nexttile()
%plot_matrix(abs(wtsg1(fsg1>fintrest,:))',tcwt,fsg1(fsg1>fintrest),'n');
surface(tcwt,fsg1(fsg1>fintrest),abs(wtsg1(fsg1>fintrest,:)));
axis tight
shading flat
colormap('jet')
ylabel('Frequency (Hz)')
nexttile()
%plot_matrix(abs(wtsg2(fsg2>fintrest,:))',tcwt,fsg2(fsg2>fintrest),'n');
surface(tcwt,fsg2(fsg2>fintrest),abs(wtsg2(fsg2>fintrest,:)));
axis tight
shading flat
ylabel('Frequency (Hz)')
colormap('jet')
xlabel('Time (ms)')

%% Modulation index analysis
%https://www.frontiersin.org/articles/10.3389/fnins.2017.00487/full#B8
%https://science.sciencemag.org/content/313/5793/1626/tab-figures-data

%low frequency phases
phase1 = angle(hilbert(oscs(1,t>100)));
phase2 = angle(hilbert(oscs(2,t>100)));

%high frequency envelopes (20-60Hz)
envelps1 = abs(wtsg1(1:281,:));
envelps2 = abs(wtsg2(1:281,:));

MI_norm11 = zeros(length(1:281),1);
MI_norm22 = zeros(length(1:281),1);
MI_norm12 = zeros(length(1:281),1);
MI_norm21 = zeros(length(1:281),1);

for ch = 1:length(MI_norm11)
    MI_norm11(ch)=abs(mean(envelps1(ch,:).*exp(1i.*phase1)));
    MI_norm12(ch)=abs(mean(envelps2(ch,:).*exp(1i.*phase1)));
    MI_norm22(ch)=abs(mean(envelps2(ch,:).*exp(1i.*phase2)));
    MI_norm21(ch)=abs(mean(envelps1(ch,:).*exp(1i.*phase2)));
end

figure()
plot(fsg1(1:281),MI_norm11,'r','LineWidth',2)
hold on
plot(fsg1(1:281),MI_norm12,'--r','LineWidth',2)
plot(fsg1(1:281),MI_norm22,'b','LineWidth',2)
plot(fsg1(1:281),MI_norm21,'--b','LineWidth',2)
xlabel('Frequency')
ylabel('Modulation Index')
legend('LF 1 on HF 1','LF 1 on HF 2','LF 2 on HF 2','LF 2 on HF 1','Location','bestoutside')
xlim([-inf,100])
%% function defitions 
% this is the 1/f noise + gaussian function for the LFP (see Materials and Methods)
function P = map(x,sig,c,amp)
    g = amp*(1/2*sqrt(pi*sig))*exp(-0.5*((x-c)/sig).^2);
    P = (1./x.^0.4).*(0.3)+g;
    %P = g;
end