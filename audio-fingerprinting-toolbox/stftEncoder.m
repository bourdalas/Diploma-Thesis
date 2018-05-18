function [keys, values] = stftEncoder(signal, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: createDatabase
% Date: 2018-05-16
% Programmer: Grigoris Bourdalas
%
% Description:
%   Function that encodes signal into audio fingerprint set of key value pairs 
%   Uses STFT for acquiring the spectrogram
%   Tiles the spectrogram, gets the local maxima of the tiles
%   Encodes them into keys and corresponding values
%   Can be used in order to encode the reference db index if store
%   mode is selected
%
% Input:
%   signal:    input audio signal 
%
% Input (optional):
%   storeFlag       flag for reference storing mode
%   storeIndex:     reference index for reference storing mode
%
% Input (optional):
%   
%
% Output:
%   keys:    encoded keys of signal
%   values:  encoded values of signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

checkPlot = @(x) any(validatestring(x, {'spectrogram','maxima','both'}));
checkMode = @(x) any(validatestring(x, {'on', 'off'}));

addRequired(p,'signal',@isnumeric);
addOptional(p,'storeFlag',checkMode)
addOptional(p,'storeIndex',0,@isnumeric)
addOptional(p,'fs',11025,@isnumeric)
addParameter(p,'plot',checkPlot)
p.KeepUnmatched = true;
parse(p,signal,varargin{:})

if strcmp(p.Results.storeFlag,'on') && (p.Results.storeIndex< 1)
    error('Storing mode Failed');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = p.Results.fs;
tfft = 0.064;
thop = 0.032;
nfft = floor(fs*tfft);
window = hamming(floor(fs*tfft));
noverlap = floor(fs*thop);

[~,F,T,P] = spectrogram(p.Results.signal,window,noverlap,nfft,fs,'yaxis');
P = 10*log(abs(P)+eps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Local Maxima
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dttile = 0.4;
dftile = 400;

DFvector = (F(length(F))-F(1))/length(F);
DTvector = (T(length(T))-T(1))/length(T);

DFindex = ceil(dftile/DFvector);
DTindex = ceil(dttile/DTvector);

PTLim = ceil(length(T)/DTindex);
PFLim = ceil(length(F)/DFindex);

tempRect = [];
rectPmax = cell(PFLim,PTLim);

DTtile = DTvector*DTindex;
DFtile = DFvector*DFindex;

maxCoords = [];
%for every tile find max point
for i=1:PFLim
    for j=1:PTLim
       if i==PFLim && j==PTLim
            tempRect = P(1+(i-1)*DFindex:size(P,1),1+(j-1)*DTindex:size(P,2));
       elseif i==PFLim
            tempRect = P(1+(i-1)*DFindex:size(P,1),1+(j-1)*DTindex:j*DTindex);
       elseif j == PTLim 
            tempRect = P(1+(i-1)*DFindex:i*DFindex,1+(j-1)*DTindex:size(P,2));
       else
            tempRect = P(1+(i-1)*DFindex:i*DFindex,1+(j-1)*DTindex:j*DTindex);
       end
       [~,I] = max(tempRect(:));
       [I_row,I_col] = ind2sub(size(tempRect),I);
       maxCoords = cat(1, maxCoords, [((j-1)*DTindex + I_col)*DTvector,((i-1)*DFindex + I_row)*DFvector]);
       rectPmax{i,j} = [((j-1)*DTindex + I_col)*DTvector,((i-1)*DFindex + I_row)*DFvector];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the Pairs of Local Maxima
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DTmax = 3;
DFmax = 350;
pair = [];
pairStartTime = [];
testP = zeros(PFLim,PTLim);

%for every rectangle 
for counter=0:PFLim*PTLim-1
    %initialize necessary temporary variables
    if counter == 0
        tempPair = [];
        tempPairStartTime = [];
    end
    tempPoints = [];

    %for the rectangles examined in time 
    for j=floor(counter/PFLim)+1:floor(counter/PFLim) + ceil(DTmax/DTtile)
       %for the rectangles examined in frequency
       for i=mod(counter,PFLim)+1:mod(counter,PFLim) +1 + ceil(DFmax/DFtile)
            %if the rectangles are in the correct space 
            if (PFLim +1 - i)> 0 && (PFLim +1 - i)<=size(rectPmax,1) && j<=size(rectPmax,2)
                testP(PFLim +1- i, j) = counter + 1; 
            end
            %collect the candidate points that may be paired
            if i<size(rectPmax,1) && j<size(rectPmax,2)
                tempPoints = cat(1,tempPoints,rectPmax{i,j});
            end
        end
    end

    %create pairs 
    if size(tempPoints,1)>1
        for tempPoint=2:length(tempPoints)
            %if the two points meet the pruning requirements
            if (tempPoints(tempPoint,1) - tempPoints(1,1))<DTmax && (tempPoints(tempPoint,1) - tempPoints(1,1))>0 && (tempPoints(tempPoint,2) - tempPoints(1,2))<DFmax
                %save the key [f1,f2,t2-t1]
                tempPair = cat(1,tempPair,[tempPoints(1,2),tempPoints(tempPoint,2),(tempPoints(tempPoint,1) - tempPoints(1,1))]);
                %save t1
                tempPairStartTime = cat(1,tempPairStartTime,tempPoints(1,1));
            end
        end
    end
end

pair = cat(1,pair,tempPair);
pairStartTime = cat(1,pairStartTime,tempPairStartTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Encode Pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%encode keys and indexes
n1 = 13;
n2 = 13;
n3 = 6;
keys = [];

ni1 = 22;
ni2 = 10;
values = [];

testdt = [];
for i=1:length(pair)
    f1 = dec2bin(floor(pair(i,1)/max(F)*(2^n1 - 1)),n1);
    f2 = dec2bin(floor(pair(i,2)/max(F)*(2^n2 - 1)),n2);
    testdt = cat(1,testdt,pair(i,3)/DTmax*(2^n3 - 1));
    dt = dec2bin(floor(pair(i,3)/DTmax*(2^n3 - 1)),n3);
    keys  = cat(1,keys,[f1 f2 dt]);
    t1 = dec2bin(floor(pairStartTime(i)/max(T)*(2^ni1 - 1)),ni1);
    if strcmp(p.Results.storeFlag,'on')
        r = dec2bin(p.Results.storeIndex,ni2);
        values = cat(1,values,[t1 r]);   
    else
        values = cat(1,values,t1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(p.Results.plot, 'spectrogram') || strcmp(p.Results.plot, 'both')
    surf(T,F,P,'edgecolor','none'); axis tight; view(0,90);grid
    colormap(autumn);
%     set(gca,'clim',[-80 -30]);
    set(gca,'xtick',[0:0.5:max(T)-eps])
    set(gca,'ytick',[0:500:max(F)-eps])
    title('Fourier Spectrogram')
    xlabel('Time (Seconds)'); ylabel('Frequency (Hz)');
    if strcmp(p.Results.plot, 'both')
        hold on;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Local Maxima
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(p.Results.plot, 'maxima') || strcmp(p.Results.plot, 'both')
    plot(maxCoords(:,1),maxCoords(:,2),'kx','MarkerSize',10,'LineWidth',1.5)
    grid on 
    ax = gca;
    ax.LineWidth = 1;
    ax.GridAlpha = 1;
    set(gca,'xtick', [0:DTtile:max(T)-2*eps]);
    set(gca,'ytick', [0:DFtile:max(F)-eps]);
    title('Fourier Spectrogram')
    xlabel('Time (Seconds)'); ylabel('Frequency (Hz)');
    set(gcf, 'renderer', 'painters')
    set(gca, 'Layer', 'top')
end


end

