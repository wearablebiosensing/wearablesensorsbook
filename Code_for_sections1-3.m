%Written in Wearable Biosensing Lab, University of Rhode Island

%Contact email: uri_wbl@ele.uri.edu

%%
%Book Chapter - Matlab Code

%sections are to but run individually

%%
%Figure 1
clear all
load('figure1data.mat')
x = data;
fs = 256;
xmax = max(abs(x));     % find the maximum value
x = x/xmax;             % scalling the signal

% time & discretisation parameters
N = length(x);
t = (0:N-1)/fs;       

% plotting of the waveform
f1 = figure('Name','Sample Data','NumberTitle','off');
f1.Color = [0.6 0.7 0.9];
f1.ToolBar = 'none';
ax1 = axes('Parent', f1);
ax1.Title.FontWeight = 'normal';
ax1.Title.FontSize = 24;
ax1.Title.FontAngle = 'italic';

plot(ax1, t, x,'red','linewidth',2)
xlim([0 max(t)])
ylim([-1.1*max(abs(x)) 1.1*max(abs(x))])
grid on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Normalized amplitude')
title('The signal in the time domain')

%%
%Section 1.1: Example 1
A = 10; 
f0 = 1000; 
phi = pi/2;
T0 = 1/f0;
tt = -2*T0 : T0/40 : 2*T0;
xx = A*cos(2*pi*f0*tt + phi);

f2 = figure('Name','Example 1','NumberTitle','off');
f2.Color = [0.6 0.7 0.9];
f2.ToolBar = 'none';
ax2 = axes('Parent', f2);
ax2.Title.FontWeight = 'normal';
ax2.Title.FontSize = 24;
ax2.Title.FontAngle = 'italic';

plot(ax2,tt,xx,'red','linewidth',2)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
title('Sinusoid: x(t) = 10 cos(2*pi*1000*t + pi/2)');
xlabel('Time (sec)');
grid on

%%
%Figure - FFT
clear all
load('figure1data.mat')
x = data;
fs = 256;
N = length(x);

xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x):fs/2;

f3 = figure('Name','Sample Data - FFT','NumberTitle','off');
f3.Color = [0.6 0.7 0.9];
f3.ToolBar = 'none';
ax3 = axes('Parent', f3);
ax3.Title.FontWeight = 'normal';
ax3.Title.FontSize = 24;
ax3.Title.FontAngle = 'italic';


plot(ax3, freq,10*log10(psdx),'red','linewidth',2)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

%%
%Section 1.2: Spectrogram
% plotting the periodogram for an even-length signal
clear all
load('figure1data.mat')
x = data;
fs = 80;

% plotting of the spectrogram

f4 = figure('Name','Sample Data - Spectrogram','NumberTitle','off');
f4.Color = [0.6 0.7 0.9];
f4.ToolBar = 'none';
ax4 = axes('Parent', f4);
ax4.Title.FontWeight = 'normal';
ax4.Title.FontSize = 24;
ax4.Title.FontAngle = 'italic';

%spectrogram(x, 1024, 3/4*1024, [], fs, 'yaxis')
spectrogram(x, [], [], [], fs, 'yaxis')
h = colorbar;
set(h, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(h, 'Magnitude, dB');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency, Hz')
title('Spectrogram of the signal')

%%
%Section 2.2.4: Example filtering - 3- point averager
clear all

nn = 0:99;			%<--Time indices
xx = cos( 0.08*pi*nn );	%<--Input signal
bb = [1/3 1/3 1/3];		%<--Filter coefficients
yy = firfilt(bb, xx);		%<--Compute the output

%%
%Section 2.2.4: Example filtering
clear all

bb = [0.5, 0.5];		%<? Filter Coefficients
ww = -pi:(pi/100):pi;	%<? omega hat

f5 = figure('Name','Filtering - Example','NumberTitle','off');
f5.Color = [0.6 0.7 0.9];
f5.ToolBar = 'none';
ax5 = axes('Parent', f5);
ax5.Title.FontWeight = 'normal';
ax5.Title.FontSize = 24;
ax5.Title.FontAngle = 'italic';
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

HH = freqz(bb, 1, ww);
subplot(2,1,1);
plot(ww, abs(HH),'red','linewidth',2)
title('2nd Order Filter - Magnitude')
subplot(2,1,2);
plot(ww, angle(HH),'red','linewidth',2)
xlabel('Normalized Radian Frequency')

title('2nd Order Filter - Linear Phase')

%%
%Section 2.2.8: Example filtering - nulling
clear all

bb = [1, -2*cos(0.5*pi), 1];		%<? Filter Coefficients
ww = -pi:(pi/100):pi;	%<? omega hat
HH = freqz(bb, 1, ww);

f6 = figure('Name','Filtering - Nulling','NumberTitle','off');
f6.Color = [0.6 0.7 0.9];
f6.ToolBar = 'none';
ax6 = axes('Parent', f6);
ax6.Title.FontWeight = 'normal';
ax6.Title.FontSize = 24;
ax6.Title.FontAngle = 'italic';
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

subplot(2,1,1);
plot(ww, abs(HH),'red','linewidth',2)
title('3rd Order Filter - Magnitude')
subplot(2,1,2);
plot(ww, angle(HH),'red','linewidth',2)
xlabel('Normalized Radian Frequency')
title('3rd Order Filter - Phase')

%%
%Section 2.3: Example filtering - simple bandpass
clear all
L = 8;
beta = 2/L; %change 2 to 1 for unity passband
fb = 64;
fs = 128;

bb = [beta*cos((2*pi*fb)/fs)]*ones(L,1);
ww = -pi:(pi/100):pi;	%<? omega hat
HH = freqz(bb, 1, ww);

f7 = figure('Name','Filtering - Simple Bandpass','NumberTitle','off');
f7.Color = [0.6 0.7 0.9];
f7.ToolBar = 'none';
ax7 = axes('Parent', f7);
ax7.Title.FontWeight = 'normal';
ax7.Title.FontSize = 24;
ax7.Title.FontAngle = 'italic';
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

subplot(2,1,1);
plot(ww, abs(HH),'red','linewidth',2)
title('3rd Order Filter - Magnitude')
subplot(2,1,2);
plot(ww, angle(HH),'red','linewidth',2)
xlabel('Normalized Radian Frequency')
title('8th Order Filter - Phase')


%%
%Section 2.3.1: Example filtering - better BPF
clear all

L = 8;
wc = 0.44*pi;
bb = zeros(L,1);
%bb = [0.54*cos(wc*((L-1)/2))  -0.46*cos((2*pi)/(L-1))*cos(wc*((L-1)/2))];		%<? Filter Coefficients
bb(1:2:end) = 0.54*cos(wc*((L-1)/2));
bb(2:2:end) = -0.46*cos((2*pi)/(L-1))*cos(wc*((L-1)/2));
ww = -pi:(pi/100):pi;	%<? omega hat
HH = freqz(bb, 1, ww);

f8 = figure('Name','Filtering - Better BPF','NumberTitle','off');
f8.Color = [0.6 0.7 0.9];
f8.ToolBar = 'none';
ax8 = axes('Parent', f8);
ax8.Title.FontWeight = 'normal';
ax8.Title.FontSize = 24;
ax8.Title.FontAngle = 'italic';
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

subplot(2,1,1);
plot(ww, abs(HH),'red','linewidth',2)
title('8th Order Filter - Magnitude')
subplot(2,1,2);
plot(ww, angle(HH),'red','linewidth',2)
xlabel('Normalized Radian Frequency')
title('8th Order Filter - Phase')



%%
%Section 3.4: Example - Symbolic
clear all

syms t A omega phi
xt = int(A*t^3 + pi, t)
yt = diff(A*cos(omega*t + phi), t)
zz = int(A*cos(omega*t + phi), t, 0, 1/omega)
zsimp = simplify(zz)


%%
%Section 3.4: Example - Voltage plot
clear all

syms v t A omega phi
v = A*cos(omega*t + phi)
v1 = subs(v, {A,omega,phi}, {100,120*pi,0})
tn = (0:0.01:1)/60;
v1n = double( subs(v1,t,tn) );  %<-- convert to numeric

f9 = figure('Name','Plotting a Symbolic Experession','NumberTitle','off');
f9.Color = [0.6 0.7 0.9];
f9.ToolBar = 'none';
ax9 = axes('Parent', f9);
ax9.Title.FontWeight = 'normal';
ax9.Title.FontSize = 24;
ax9.Title.FontAngle = 'italic';
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

plot(ax9, tn, v1n,'red','linewidth',2)
ylabel('Voltage')
xlabel('Time')
title('Voltage Plot')

%%
%Section 3.5: Example Fourier
clear all
syms wt t ak k
N=10; T0=2;
ak = sin(k)/k;
wt = fouriersynth(ak, N, T0);

f10 = figure('Name','Example of Fourier Series','NumberTitle','off');
f10.Color = [0.6 0.7 0.9];
f10.ToolBar = 'none';
ax10 = axes('Parent', f10);
ax10.Title.FontWeight = 'normal';
ax10.Title.FontSize = 24;
ax10.Title.FontAngle = 'italic';
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

fplot(wt, [-T0,T0], 'red', 'linewidth', 2); grid on
ylabel('Voltage')
xlabel('Time')
axis tight 			%<---make fplot show the whole thing
title('Example of Fourier Series')
