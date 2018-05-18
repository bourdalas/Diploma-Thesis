function [keys, values] = cqtEncoder(signal, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: createDatabase
% Date: 2018-05-16
% Programmer: Grigoris Bourdalas
%
% Description:
%   Function that encodes signal into audio fingerprint set of key value pairs 
%   Uses CQT for acquiring the spectrogram
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
addParameter(p,'plot','off',checkPlot)
p.KeepUnmatched = true;
parse(p,signal,varargin{:})

if strcmp(p.Results.storeFlag,'on') && (p.Results.storeIndex< 1)
    error('Storing mode Failed');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = p.Results.fs;
fmin = 51.91;
r = 3;
B = r*12;
gamma = 0; 
fmax = fs/2;
signalLen = length(signal);

Xcq = cqt(signal, B, fs, fmin, fmax, 'rasterize', 'full','gamma', gamma, 'winfun', 'hamming', 'phasemode', 'local');
C = 10*log10(abs(Xcq.c)+eps);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Local Maxima
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = linspace(0,signalLen/fs,size(C,2));
F = [1:size(C,1)]';

Fy = Xcq.fbas;

dttile = 0.4;
DBins = 18;

DTvector = (T(length(T))-T(1))/length(T);
DTindex = ceil(dttile/DTvector);

PTLim = floor(length(T)/DTindex);
PFLim = length(1:DBins:length(F));

rectPmax = cell(PFLim,PTLim);

DTtile = DTvector*DTindex;
maxCoords = [];
%for every tile find max point
for i=1:PFLim
    for j=1:PTLim
       if i*DBins-1 > length(F)
           tempRect = C((i-1)*DBins+1:length(F),1+(j-1)*DTindex:j*DTindex);
       else
           tempRect = C((i-1)*DBins+1:i*DBins-1,1+(j-1)*DTindex:j*DTindex);
       end
       [~,I] = max(tempRect(:));
       [I_row,I_col] = ind2sub(size(tempRect),I);
       maxcordx = ((j-1)*DTindex + I_col -1)*DTvector;
       if (i-1)*DBins + I_row > size(C,1)
           maxcordy = F(end);
       else
           maxcordy = F((i-1)*DBins + I_row);
       end
       maxCoords = cat(1, maxCoords, [maxcordx Fy(maxcordy)]);
       rectPmax{i,j} = [maxcordx,maxcordy];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the Pairs of Local Maxima
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DTmax = 1.2;
DBmax = 24;
pair = [];
pairStartTime = [];
testP = zeros(PFLim,PTLim);
tic;
%for every rectangle 
for counter=0:PFLim*PTLim-1
    %initialize necessary temporary variables
    if counter == 0
        tempPair = [];
        tempPairStartTime = [];
    end
    tempPoints = [];

    %for the rectangles examined in time 
    for j=floor(counter/PFLim)+1:floor(counter/PFLim)+1 + ceil(DTmax/DTtile)
       %for the rectangles examined in frequency
       for i=mod(counter,PFLim)+1:mod(counter,PFLim)+1 + DBmax
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
            if (tempPoints(tempPoint,1) - tempPoints(1,1))<DTmax && (tempPoints(tempPoint,1) - tempPoints(1,1))>=0 && (tempPoints(tempPoint,2) - tempPoints(1,2))<DBmax && (tempPoints(tempPoint,2) - tempPoints(1,2))>=0
                %save the key [floor(b1/6),b2-b1,t2-t1]
                f1 = F(tempPoints(1,2));
                f2 = F(tempPoints(tempPoint,2));
                tempPair = cat(1,tempPair,[floor(f1/6),floor(f2 - f1),(tempPoints(tempPoint,1) - tempPoints(1,1))]);
%                     tempPair = cat(1,tempPair,[floor(tempPoints(1,2)/6),tempPoints(tempPoint,2) - tempPoints(1,2),(tempPoints(tempPoint,1) - tempPoints(1,1))]);
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
    surf(T,F,C,'edgecolor','none'); axis tight; view(0,90);grid
    xlabel('Time (seconds)', 'FontSize', 12, 'Interpreter','latex'); 
    ylabel('Frequency (Herz)', 'FontSize', 12, 'Interpreter','latex');
    set(gca,'xtick', [0:DTtile:max(T)])
    set(gca,'ytick', Fy(1:DBins:end))
    set(gca,'clim',[min(C(:)); max(C(:));]);
    title('CQT Spectrogram')
    colormap(autumn)
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
    set(gca,'xtick',[0:DTtile:max(T)])
    set(gca,'ytick',Fy(1:DBins:end))
    title('CQT Spectrogram')
    xlabel('Time (seconds)', 'FontSize', 12, 'Interpreter','latex'); 
    ylabel('Frequency (Herz)', 'FontSize', 12, 'Interpreter','latex');
    set(gcf,'renderer', 'painters')
    set(gca,'Layer', 'top')
end

