clear all;
clc

createDbFlag = 0;
identifyRandomExcerptFlag = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Database Creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if createDbFlag
    rootPath = '/media/greg/BE06B21E06B1D81B/Users/Greg4/Music/complete/';
    dbName = 'STFT_50';
    dbSize = 50;
    method = 'STFT';
    createDatabase(rootPath, dbSize, 'verbose', 'on','dbName', dbName, 'method', method);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Random Excerpt Identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if identifyRandomExcerptFlag
    % Specify dbName if script used only for identification
    if ~createDbFlag  
        dbName = 'STFT_50';
        method = 'STFT';
    end
    load(fullfile(pwd,'Databases',dbName,strcat(dbName,'.mat')));
    randID = floor((dbInfo.dbSize-1)*rand(1,1)+1);
    excerptFilePath = fullfile(dbInfo.paths{randID},dbInfo.names{randID});
    audioStart = 10;
    audioLength = 30;
    P = 6;
    subLength = 5;

    disp(['Searching for:   ', dbInfo.names{randID}])
    disp(['Starting at: ', num2str(audioStart), ' seconds'])

    [excerpt, fidelityFlag] = audioFilePreprocessor(excerptFilePath,'audioStart',audioStart,'audioLength',audioLength);

    if fidelityFlag
        [referenceIndex, referenceStartinTime] = identifyExcerpt(excerpt, dbName, P, subLength, 'method', method);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display Results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    verbose_mode_result = 'on';
    if strcmp(verbose_mode_result, 'on') 
        if ~isequal(referenceIndex,-1)
            disp('.');
            disp(['Reference Found: ' dbInfo.names{referenceIndex,1}]);
            disp(['Starting at: ' num2str(referenceStartinTime) ' seconds']);
        else
            disp('The excerpt was not found in the database');
        end
        disp('.');
        disp('.');
    end
end