function [audioSignal, fidelityFlag] = audioFilePreprocessor(filename, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: createDatabase
% Date: 2018-05-16
% Programmer: Grigoris Bourdalas
%
% Description:
%   Function that preprocesses an audio file. 
%   Adjusts the sample rate in to a specific one set in the input args
%   Converts multiChannel or stereo audio files to mono
%   Crops the signal as requested 
%
% Input:
%   filename:    path of the audio file
%
% Input (optional):
%   audioStart:     the starting second of the references
%   audioLength:    the length in seconds of the references stored
%   fs:             desired sampling frequency
%
% Output:
%   audioSignal:    the output audio signal
%   fidelityFlag:   true if file is pre-processed succesfully, false otherwise 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
validMethods = {'WANG','CQT'};
checkMethod = @(x) any(validatestring(x, validMethods));

addRequired(p,'filename',@isfile);
addOptional(p,'audioStart',@isnumeric)
addOptional(p,'audioLength',@isnumeric)
addOptional(p,'fs',11025,@isnumeric)
p.KeepUnmatched = true;
parse(p,filename,varargin{:})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre Process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[audioSignal,initFs] = audioread(p.Results.filename);

%Adjust Sampling Frequency
audioSignal = audioSignal(1:floor(initFs/p.Results.fs):end,:);
%Convert Stereo to Mono
audioSignal = sum(audioSignal, 2) / size(audioSignal, 2);

%Crop the signal if needed and possible
duration = floor(length(audioSignal)/p.Results.fs);
fidelityFlag = 0;
sampleStart = p.Results.audioStart*p.Results.fs + 1;
sampleEnd = p.Results.audioStart*p.Results.fs + p.Results.audioLength*p.Results.fs;
if (duration > p.Results.audioStart + p.Results.audioLength)
    audioSignal = audioSignal(sampleStart:sampleEnd);
    fidelityFlag = 1;
end

end

