clear all;

%% PARAMETERS
fs = 11025;
fmin = 51.9;
B = 36;
gamma = 0; 
fmax = fs/2;

shiftBins = 36;

%% INPUT SIGNAL
x = audioread('kempff1.wav');
%x = wavread('kempff1.wav');
x = x(:); xlen = length(x);

%% COMPUTE COEFFIENTS
Xcq = cqt(x, B, fs, fmin, fmax, 'rasterize', 'full', 'gamma', gamma, ...
    'phasemode', 'local', 'normalize' , 'sine', 'winfun', 'hann');
c = Xcq.c;
figure(1)
plotnsgtf(Xcq.c,Xcq.shift, fs, fmin, fmax, B, 2, 120)
% %% PITCH SHIFTING
% if shiftBins ~= 0
%     Y = phaseUpdate(c,Xcq.fbas,shiftBins,Xcq.xlen, fs, 1e-6);
%     Y = circshift(Y,shiftBins);
%     if shiftBins > 0
%         Y(1:shiftBins) = 0;
%     elseif shiftBins < 0
%         Y(end-shiftBins+1:end) = 0;
%     end
%     Xcq.c = Y;
% end
% 
% %% ICQT
% [y gd] = icqt(Xcq);
% 
% plot(x)
% hold on 
% plot(y,'r')
% hold on
% plot(real(c))