%% 
%This file runs the code to generate the simulations in Fig. S1. It plots
%the LFPs phase shifted as well as the FFT spectrum
dt = 0.05;
t = 0:dt:2000;
fs = 1000/dt;
L = length(t);
f = fs*(0:(L/2))/L;

fdiff = diff(f);
fstep = fdiff(1);

signal1 = zeros(1,length(t));
signal2 = zeros(1,length(t));

fLow = 0.1;
fHigh = 1000;
for ff = fLow:fstep:fHigh 
    phase1 = pi*(2*rand()-1);
    phase2 = phase1 - (pi/2+0.3*pi*(2*rand()-1));%this is the phase difference DELTA THETA (see materials and methods)
    %signal = signal + (1/ff)*sin(2*pi*(ff/1000)*t-phase);
    signal1 = signal1 + map(ff,0.5,8,10)*sin(2*pi*(ff/1000)*t-phase1);
    signal2 = signal2 + map(ff,0.5,8,10)*sin(2*pi*(ff/1000)*t-phase2);
end

P1 = fftshift(fft(signal1));
P1 = P1.*conj(P1);
P2 = fftshift(fft(signal2));
P2 = P2.*conj(P2);

tiledlayout(2,1)
nexttile()
%time series plot
plot(t,signal1,'k','LineWidth',2)
hold on
plot(t,signal2,'r','LineWidth',2)
hold off
xlabel('Time (ms)')
ylabel('Amplitude (a.u.)')
nexttile()
%FFT plot, log-scale to show 1/f noise
plot(log10(f),10*log10(abs(P1(ceil(length(P1)/2):end))),'k','LineWidth',2)
hold on 
plot(log10(f),10*log10(abs(P2(ceil(length(P2)/2):end))),'r','LineWidth',2)
hold off
xlim([0 log10(1000)])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
xticks([0 1 2 3])
xticklabels({'0','10','100','1000'})
%plot(f,P(ceil(length(P)/2):end))
%xlim([0 100])

%% this is the 1/f noise + gaussian function for the LFP (see Materials and Methods)
function P = map(x,sig,c,amp)
    g = amp*(1/2*sqrt(pi*sig))*exp(-0.5*((x-c)/sig).^2);
    P = (1./x.^0.4).*(0.3)+g;
end