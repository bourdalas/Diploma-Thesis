function [finalReferenceIndex, finalStartingTime] = identifyExcerpt(excerpt, dbName, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: createDatabase
% Date: 2018-05-16
% Programmer: Grigoris Bourdalas
%
% Description:
%   Function that identifies an unknown excerpt by searching in the input
%   database. it splits the signal into P subsignals of length set by the 
%   user. For each subsignal, the encoding takes place, its keys are 
%   searched in the database, the best match is found. The final ouput
%   is decided by examing the best matches with the fusion of local
%   decisions module.
%
% Input:
%   signal:    input audio signal 
%
% Input (optional):
%   P:                   number of splited sub-signals 
%   subsignalLength:     length of the subsignals
%
% Output:
%   keys:    encoded keys of signal
%   values:  encoded values of signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

checkPlot = @(x) any(validatestring(x, {'spectrogram','maxima','both'}));
checkMethod = @(x) any(validatestring(x, {'STFT','CQT'}));


addRequired(p,'excerpt',@isnumeric);
addRequired(p,'dbName',@ischar);
addOptional(p,'P',1,@isnumeric)
addOptional(p,'subsignalLength',5,@isnumeric)
addParameter(p,'plot',checkPlot)
addParameter(p,'fs',11025,@isnumeric)
addParameter(p,'method','CQT',checkMethod)
p.KeepUnmatched = true;
parse(p,excerpt,dbName, varargin{:})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excerpt Splitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = p.Results.P;
subsignalLength = p.Results.subsignalLength;
fs = p.Results.fs;
overlap = 0.5;

subsignals = cell(P,1);
for ii=1:P
    subsignals{ii} = excerpt(1 + floor(fs*(ii-1)*ceil(overlap*subsignalLength)):ceil(fs*((ii-1)*ceil(overlap*subsignalLength) + subsignalLength)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(pwd,'Databases',dbName,strcat(dbName,'.mat')));
dbInfo.referenceLength = 60;
database = lmdb.DB(strcat('./', fullfile('Databases', dbName)), 'MAPSIZE', 1024^3, 'RDONLY', true, 'DUPSORT', true,'NOLOCK', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Match Subsignals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = zeros(P,1);
s = zeros(P,1);
Tvote = ceil(P/2);
BinLength = 1;
ni1 = 22;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each subsignal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cntr=1:P
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Encode Subsignal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    signal = subsignals{cntr};
    if strcmp(p.Results.method,'CQT')
        [excerptkeys, excerptt1] = cqtEncoder(signal,'storeFlag', 'off');
    elseif strcmp(p.Results.method, 'STFT')
        [excerptkeys, excerptt1] = stftEncoder(signal, 'storeFlag', 'off');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Search for matching keys 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    histograms = zeros(dbInfo.dbSize,dbInfo.referenceLength/BinLength);
    %decode starting time of keys back in decimals
    excerptt1 = bin2dec(excerptt1)*subsignalLength/(2^ni1-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each key in the excerpt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for i=1:size(excerptkeys,1)
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fetch its values from the database
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        cursor = database.cursor('RDONLY', true);
        key = [];
        values = [];

        if cursor.find(excerptkeys(i,:))
            key = cursor.key;
            values = cursor.value;
        end

        while length(key)==32 && cursor.next()  
            key = cursor.key;
            if strcmp(key, excerptkeys(i,:))    
                values = cat(1,values, cursor.value);
            else
                key = [];
                break;
            end
        end
        clear cursor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update Histogram array
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        for valueCounter=1:size(values,1)
            if size(values,2) == 32
                % Decode starting time of reference key in decimals
                referenceTime = bin2dec(values(valueCounter,1:ni1))*dbInfo.referenceLength/(2^ni1-1);
                % Do the incremental operation in the reference's histogram
                if referenceTime - excerptt1(i) > 0
                    referenceIndex = bin2dec(values(valueCounter,ni1+1:32));
                    d = round(referenceTime - excerptt1(i));
                    histograms(referenceIndex,d+1) = histograms(referenceIndex,d+1) + 1;
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %Search for the maximum peak in the histogram array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    I = zeros(dbInfo.dbSize,1);
    histPeak = zeros(dbInfo.dbSize,1);
    for i=1:dbInfo.dbSize
        [histPeak(i),I(i)] = max(histograms(i,:));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detect reference index and the excerpt's starting time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    [~,referenceIndex] = max(histPeak);
    s(cntr) = I(referenceIndex)-1;    
    r(cntr) = referenceIndex;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fusion of the local decisions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EVAL = zeros(P,P);
finalStartingTimeArray = zeros(P,1);
finalReferenceIndex = -1;
finalStartingTime = -1;

for i=1:P
    finalStartingTimeArray(i) = ceil(s(i) - (i-1)*ceil(subsignalLength*(1 - overlap)));
    for j=1:P
        if j ~= i
             EVAL(i,j) = floor(s(i) - (i-1)*ceil(subsignalLength*(1 - overlap))) == floor(s(j) - (j-1)*ceil(subsignalLength*(1 - overlap))) && r(i) == r(j); 
        else
            EVAL(i,j) = 1;
        end            
    end
end

for i=1:P
   if size(find(EVAL(:,i) == 1),1) >= Tvote 
      finalReferenceIndex = r(i);
      finalStartingTime = finalStartingTimeArray(i); 
      break;
   end
end

end