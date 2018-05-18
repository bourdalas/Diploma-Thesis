function [f_audio_out,timepositions_afterDegr] = degradationUnit_adaptiveEqualizer(f_audio, samplingFreq, timepositions_beforeDegr, parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: degradationUnit_adaptiveEqualizer
% Date of Revision: 2013-01-23
% Programmer: Sebastian Ewert
%
% Description:
% - designs a filter such that the mean spectrum of f_audio becomes similar
% to a given mean spectrum
% - mean spectra are specified using a decibel scale, i.e. if x is a magnitude
% spectrum, then apply x -> 20*log10(x)
% - there are four ways to specify the destination mean spectrum: 1. by
% loading example files provided with the toolbox, 2. by using specific
% noise "color" profiles, 3. by providing the destination mean spectrum
% using the parameter destMagFreqResp, 4. by providing audio data from
% which the destination mean spectrum is computed.
% - if mean spectra are computed, this is done by computing a magnitude
% spectrogram, deriving the corresponding values in dB, and then averaging
% over all time frames.
% - the same is done for f_audio and the two mean spectral vectors are
% compared
% - Then a filter is designed such that f_audio's mean spectral vector
% becomes similar to destMagFreqResp
% - filtering is done using a linear-phase FIR filter 
% - this unit normalizes the output
%
% Notes:
% - note that the mean spectrum is often a good approximation of what is
% sometimes referred to as the mean spectral shape (however, this term is not
% defined for general sound recordings). In this sense, the function
% modifies the mean spectral shape to the desired one.
%
% Input:
%   f_audio      - audio signal \in [-1,1]^{NxC} with C being the number of
%                  channels
%   samplingFreq - sampling frequency of f_audio
%   timepositions_beforeDegr - some degradations delay the input signal. If
%                             some points in time are given via this
%                             parameter, timepositions_afterDegr will
%                             return the corresponding positions in the
%                             output. Set to [] if unavailable. Set f_audio
%                             and samplingFreq to [] to compute only
%                             timepositions_afterDegr.
%
% Input (optional): parameter
%   .loadInternalMagFreqResp=1  - loads one of the destMagFreqResp provided
%                                 by the toolbox
%   .internalMagFreqResp='Beatles_NorwegianWood'
%   .computeMagFreqRespFromAudio - computes destMagFreqResp from given
%                                  audio data
%   .computeMagFreqRespFromAudio_audioData
%                               - audio data for .computeMagFreqRespFromAudio
%   .computeMagFreqRespFromAudio_sf - sampl freq for .computeMagFreqRespFromAudio
%   .destMagFreqResp = []       - in db. See above.
%   .destMagFreqResp_freqs = [] - must have same length as destMagFreqResp.
%                                 In Hertz
%   .fftLength = 2 ^ nextpow2(0.02 * samplingFreq); - fft length to
%                                 calculate spectrogram of f_audio.
%
% Output:
%   f_audio_out   - audio signal \in [-1,1]^{NxC} with C being the number
%                   of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Audio Degradation Toolbox
%
% Centre for Digital Music, Queen Mary University of London.
% This file copyright 2013 Sebastian Ewert, Matthias Mauch and QMUL.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.  See the file
% COPYING included with this distribution for more information.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
    parameter=[];
end
if nargin<3
    timepositions_beforeDegr=[];
end
if nargin<2
    error('Please specify input data');
end

if isfield(parameter,'loadInternalMagFreqResp')==0
    parameter.loadInternalMagFreqResp = 1;
end
if isfield(parameter,'loadNoiseColorPreset')==0
    parameter.loadNoiseColorPreset = 0;
end
if isfield(parameter,'computeMagFreqRespFromAudio')==0
    parameter.computeMagFreqRespFromAudio = 0;
end

if isfield(parameter,'internalMagFreqResp')==0
    parameter.internalMagFreqResp = 'Beethoven_Appasionata_Rwc';
end
if isfield(parameter,'noiseColorPreset')==0
    parameter.noiseColorPreset = 'pink';
end
if isfield(parameter,'computeMagFreqRespFromAudio_audioData')==0
    parameter.computeMagFreqRespFromAudio_audioData = [];
end
if isfield(parameter,'computeMagFreqRespFromAudio_sf')==0
    parameter.computeMagFreqRespFromAudio_sf = [];
end
if isfield(parameter,'destMagFreqResp')==0
    parameter.destMagFreqResp = [];
end
if isfield(parameter,'destMagFreqResp_freqs')==0
    parameter.destMagFreqResp_freqs = [];
end
if isfield(parameter,'fftLength')==0
    parameter.fftLength = 2 ^ nextpow2(0.02 * samplingFreq);
end
if isfield(parameter,'filterOrder')==0
    parameter.filterOrder = round(parameter.fftLength/2);
end
if isfield(parameter,'visualizations')==0
    parameter.visualizations = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(f_audio)
    % we design a linear-phase filter. Such filters have the property that
    % the group delay for every frequency is always equal to
    % parameter.filterOrder/2. Therefore, although we built a signal
    % adaptive filter, we can adjust for that delay. That means that the
    % timepositions_beforeDegr don't require any adjustment.
    timepositions_afterDegr = timepositions_beforeDegr;
    return;
end

% load/compute destMagFreqResp
[destMagFreqResp,destMagFreqResp_freqs] = internal_getDestMagFreqResp(parameter);

% compute mean spectral vector for f_audio
[meanMagSpec,meanMagSpec_freqs] = internal_computeMeanSpectralVector(f_audio,samplingFreq,parameter.fftLength);
meanMagSpec = internal_standardizeMagFreqResp(meanMagSpec);

% compute magnitude response for the filter to be designed
destMagFreqResp_org = destMagFreqResp;
destMagFreqResp_freqs_org = destMagFreqResp_freqs;
if ~((length(destMagFreqResp_freqs)==length(meanMagSpec_freqs)) && (all(meanMagSpec_freqs(:) == destMagFreqResp_freqs(:))))
    % in this case we interpolate the frequency response using a
    % spline interpolation
    destMagFreqResp = spline(destMagFreqResp,destMagFreqResp_freqs,meanMagSpec_freqs);
end
filter_magFreqResp = destMagFreqResp(:) - meanMagSpec(:);

% design filter (fir2 is linear phase)
filter_magFreqResp_linear = 10 .^ (filter_magFreqResp/20);
b = fir2(parameter.filterOrder,meanMagSpec_freqs/(samplingFreq/2),filter_magFreqResp_linear);

% apply filter
parameterApplyImpulseResponse.loadInternalIR = 0;
parameterApplyImpulseResponse.impulseResponse = b;
parameterApplyImpulseResponse.impulseResponseSampFreq = samplingFreq;
parameterApplyImpulseResponse.normalizeOutputAudio = 1;
parameterApplyImpulseResponse.averageGroupDelayOfFilter = round(parameter.filterOrder/2);
[f_audio_out,timepositions_afterDegr] = degradationUnit_applyImpulseResponse(f_audio, samplingFreq, timepositions_beforeDegr, parameterApplyImpulseResponse);

if parameter.visualizations
    fvtool(b,1);
    
    [meanMagSpecOut,meanMagSpecOut_freqs] = internal_computeMeanSpectralVector(f_audio_out,samplingFreq,parameter.fftLength);
    meanMagSpecOut = internal_standardizeMagFreqResp(meanMagSpecOut);
    
    figure;
    plot(destMagFreqResp_freqs_org,destMagFreqResp_org,'y');
    hold on;
    plot(meanMagSpecOut_freqs,meanMagSpecOut,'k');
    title('Comparison: destMagFreqResp(y) and mean spectral vector of output(k)')
end


end

function [f_meanmagspec_db,freqs] = internal_computeMeanSpectralVector(f_audio,fs,fftLength)

f_audio = mean(f_audio,2);

[f_spec,freqs,time] = spectrogram(f_audio,hanning(fftLength),fftLength/2,fftLength,fs);

f_magspec_db = 20 * log10(abs(f_spec));

f_magspec_db(:,isinf(sum(abs(f_magspec_db),1))) = []; % ignore columns with -inf/inf entries
f_magspec_db(:,isnan(sum(abs(f_magspec_db),1))) = [];

f_meanmagspec_db = mean(f_magspec_db,2);

end

function magFreqResp = internal_standardizeMagFreqResp(magFreqResp)

temp = magFreqResp(~isinf(magFreqResp));
temp = temp(~isnan(magFreqResp));
maxRobust = max(temp);

magFreqResp = magFreqResp - maxRobust;

magFreqResp(magFreqResp > 0) = 0;  % remaining positive inf
magFreqResp(magFreqResp < -80) = -80;  % remaining positive inf

end

function [destMagFreqResp,destMagFreqResp_freqs,fftLength] = internal_getDestMagFreqResp(parameter)
if parameter.loadInternalMagFreqResp
    % load example included in toolbox
    
    fullFilenameMfile = mfilename('fullpath');
    [pathstr,name,ext] = fileparts(fullFilenameMfile);
    dirRootIRs = fullfile(pathstr,'../degradationData');
    
    names_internal = {'Beatles_NorwegianWood','Beethoven_Appasionata_Rwc'};
    indexInternal = find(strcmpi(names_internal,parameter.internalMagFreqResp), 1);
    if isempty(indexInternal)
        error('Please specify a valid internal name')
    end
    
    switch indexInternal
        case 1
            file = fullfile(dirRootIRs,'SpecEnvelopes/Beatles_NorwegianWood.mat');
        case 2
            file = fullfile(dirRootIRs,'SpecEnvelopes/Beethoven_Appasionata_Rwc.mat');
    end
    load(file, 'destMagFreqResp', 'destMagFreqResp_freqs');
    
elseif parameter.loadNoiseColorPreset
    switch(lower( parameter.noiseColorPreset))
        case 'white'
            freqExponent = 0;
        case 'pink'
            freqExponent = 0.5;
        case 'brown'
            freqExponent = 1;
        case 'blue'
            freqExponent = -0.5;
        case 'violet'
            freqExponent = -1;
    end
    
    lengthMagResp = parameter.fftLength/2+1;
    destMagFreqResp_freqs = linspace(0,samplingFreq/2,lengthMagResp);
    magResp = 1./destMagFreqResp_freqs.^freqExponent;
    magResp(1) = 1;
    destMagFreqResp = 20 * log10(magResp);
    
elseif parameter.computeMagFreqRespFromAudio
    % compute destMagFreqResp as mean spectral vector from given audio data
    [destMagFreqResp,destMagFreqResp_freqs] = internal_computeMeanSpectralVector(...
        parameter.computeMagFreqRespFromAudio_audioData,parameter.computeMagFreqRespFromAudio_sf,parameter.fftLength);
else
    destMagFreqResp = parameter.destMagFreqResp;
    destMagFreqResp_freqs = parameter.destMagFreqResp_freqs;
end

% standardize destMagFreqResp
destMagFreqResp = internal_standardizeMagFreqResp(destMagFreqResp);
end






