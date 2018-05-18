function createDatabase(root, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: createDatabase
% Date: 2018-05-16
% Programmer: Grigoris Bourdalas
%
% Description:
%   Function that creates a LMDB database that contains all wav files
%   inside the folders of a root directory. The wav files are preprocessed
%   accordingly, encoded and stored in the database. A database info file
%   is also created.
%
% Input:
%   root:    path of the reference files' root folder
%
% Input (optional):
%   dbSize:         preconfigure database size 
%   audioStart:     the starting second of the references
%   audioLength:    the length in seconds of the references stored
%   method:         the database method name, "STFT" or "CQT"
%   verbose:        verbose mode flag, takes the values: 'on','off'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;

checkMethod = @(x) any(validatestring(x, {'STFT','CQT'}));
checkVerbose = @(x) any(validatestring(x, {'on','off'}));

addRequired(p,'root',@isfolder);
addOptional(p,'dbSize',@isnumeric)
addOptional(p,'audioStart',0,@isnumeric)
addOptional(p,'audioLength',60,@isnumeric)
addOptional(p,'fs',11025,@isnumeric)
addParameter(p,'verbose','off',checkVerbose)
addParameter(p,'dbName',@ischar)
addParameter(p,'method','CQT',checkMethod)
parse(p,root,varargin{:})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract Reference Titles List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(root));
albumeDirs = dir(p.Results.root);
albumeDirs(~[albumeDirs.isdir]) = [];  %remove non-directories

albumDirsInfo = cell(length(albumeDirs),1);
for ii = 1 : length(root)
  albumDirsInfo{ii} = dir(cat(2,albumeDirs(ii).folder,'/',albumeDirs(ii).name,'/*.wav'));
end

refFileList = [];
for ii = 1 :length(albumDirsInfo)
 refFileList = cat(1,refFileList, albumDirsInfo{ii}(:,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference Title Files PreProcessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(p.Results.verbose, 'on')
    disp('Pre Processing reference files')
end

dbInsertList = {};
for ii=1:length(refFileList)
    refFilePath = fullfile(refFileList(ii,1).folder, refFileList(ii,1).name);
    [refSignal, fidelityFlag] = audioFilePreprocessor(refFilePath,p.Results.audioStart, p.Results.audioLength, p.Results.fs);
    if fidelityFlag
        dbInsertList = cat(1,dbInsertList,{refSignal, refFileList(ii,1).name, refFileList(ii,1).folder}); 
    end
    if (ii >= p.Results.dbSize)
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create Database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dbName = strcat(p.Results.method, '_', num2str(p.Results.dbSize));
dbSize = size(dbInsertList,1);
db = lmdb.DB(fullfile(pwd,'Databases', dbName), 'MAPSIZE', 1024^3, 'DUPSORT', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference Storing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:dbSize

    if strcmp(p.Results.verbose, 'on')
        disp(strcat('Processing reference ID: ', num2str(ii)))
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Encode Signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(p.Results.method,'CQT')
        [keys,values] = cqtEncoder(dbInsertList{ii,1}, 'storeFlag', 'on', 'storeIndex', ii);
    elseif strcmp(p.Results.method, 'STFT')
        [keys,values] = stftEncoder(dbInsertList{ii,1}, 'storeFlag', 'on', 'storeIndex', ii);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Store in database
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    transaction = db.begin();
    try
        for jj=1:size(values,1)
            transaction.put(keys(jj,:),values(jj,:));
        end
        transaction.commit();
    catch exception
        transaction.abort();
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Store database information file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dbInfo.names = dbInsertList(:,2);
dbInfo.paths = dbInsertList(:,3);
dbInfo.dbSize = dbSize;
dbInfo.referenceLength = p.Results.audioLength;

save(strcat(dbName,'.mat'), 'dbInfo');
movefile(strcat(dbName,'.mat'), fullfile('Databases', dbName));
end

