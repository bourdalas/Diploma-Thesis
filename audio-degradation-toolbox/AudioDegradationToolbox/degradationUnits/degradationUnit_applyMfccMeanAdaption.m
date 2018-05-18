function [f_audio_out,timepositions_afterDegr] = degradationUnit_applyMfccMeanAdaption(f_audio, samplingFreq, timepositions_beforeDegr, parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: degradationUnit_applyMfccMeanAdaption
% Date of Revision: 2013-01-23
% Programmer: Sebastian Ewert
%
% Description:
% - designs a filter such that the mean MFCCs of f_audio become similar
% to a given mean MFCC
% - the destination MFCC can either be provided or computed from other
% audio data.
% - since there is not just one way to compute MFCCs, the provided MFCCs
% should be computed exactly as done in the code below
% - filtering is done using a linear-phase FIR filter
% - this unit normalizes the output
%
% Notes:
% - the mean MFCC corresponds to the maximum likelihood estimate of the mean
% vector for a multivariate Gaussian fitted to the set of MFCCs
% - that means that using this function the recording in f_audio is altered
% such that it looks almost identical to the one in audioDataForDestinationMfcc
% (or the one corresponding to the mean MFCC) in terms of the mean of a
% fitted Gaussian.
% - This is intersting in the context of studying the behaviour of methods
% that employ statistics over MFCCs to classify/analyse recordings.
% - The method essentially aims for having an MFCC implementation very
% similar to the one in Malcolm Slaney's auditory toolbox.
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

if isfield(parameter,'destinationMeanMfcc')==0
    parameter.destinationMeanMfcc = [];
    % make sure to provide the same coefficients as specified in
    % parameter.MFCC_cepstralCoeffsToKeep. Also use the same method to
    % compute the MFCCs as used here.
end
if isfield(parameter,'audioDataForDestinationMfcc')==0
    parameter.audioDataForDestinationMfcc = [];
end
if isfield(parameter,'audioDataForDestMfcc_sf')==0
    parameter.audioDataForDestMfcc_sf = [];
end
if isfield(parameter,'TikhonovLambdaCandidates')==0
    parameter.TikhonovLambdaCandidates = 10.^-[10:-1:0];
end
if isfield(parameter,'TikhonovMaxError')==0
    parameter.TikhonovMaxError = 10^-5;
end
if isfield(parameter,'filterOrder')==0
    parameter.filterOrder = 500;
end
if isfield(parameter,'visualizations')==0
    parameter.visualizations = 0;
end


% All following parameters were taken from the auditory toolbox. Some of
% these value only make sense for speech but are left as is as they are
% used like this in many other toolboxes
if isfield(parameter,'MFCC_lowestFrequencyHz')==0
    parameter.MFCC_lowestFrequencyHz = 133.3333;
    % this is not a good choice for music. Bass is filtered out completely.
    % This chould be lowered to 66.66666666 and MFCC_numLinearFilters
    % increased by one.
end
if isfield(parameter,'MFCC_numLinearFilters')==0
    parameter.MFCC_numLinearFilters = 13;
end
if isfield(parameter,'MFCC_numLogFilters')==0
    parameter.MFCC_numLogFilters = 27;
end
if isfield(parameter,'MFCC_linearFreqSpacingHz')==0
    parameter.MFCC_linearFreqSpacingHz = 66.66666666;
end
if isfield(parameter,'MFCC_logFreqSpacing')==0
    parameter.MFCC_logFreqSpacing = 1.0711703;
end
if isfield(parameter,'MFCC_windowSizeSec')==0
    parameter.MFCC_windowSizeSec = 0.016;  % AT sets 256 samples for a 16kHz signal
end
if isfield(parameter,'MFCC_featuresPerSecond')==0
    parameter.MFCC_featuresPerSecond = 100;
end
if isfield(parameter,'MFCC_applyPreemphasis')==0
    parameter.MFCC_applyPreemphasis = 1;
end
if isfield(parameter,'MFCC_preemphasisFilter')==0
    parameter.MFCC_preemphasisFilter = [1 -0.97];
    % this was never adapted to work for every sampling rate! It is too
    % weak for sampling rates higher than 16000.
end
if isfield(parameter,'MFCC_cepstralCoeffsToKeep')==0
    parameter.MFCC_cepstralCoeffsToKeep = [2:13];
    % in contrast to the AT we ignore the first MFCC coefficient (DC) by default
end
if isfield(parameter,'MFCC_additionalFilter')==0
    parameter.MFCC_additionalFilter = [];
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

% compute average MFCC which we will try to achieve for f_audio as well.
% [f_spec,freqinfo,timeinfo] = spectrogram(mean(parameter.audioDataForDestMfcc,2),hanning(parameter.fftLength),parameter.fftLength/2,parameter.fftLength,parameter.audioDataForDestMfcc_sf);
[f_mfccDest,sideinfo] = internal_computeMfccs(parameter.audioDataForDestMfcc,parameter.audioDataForDestMfcc_sf,parameter);
f_mfccDestMean = internal_infNanRobustMatrixMean(f_mfccDest,2);

% compute log-mel-spectra of input
f_mfccMean = internal_infNanRobustMatrixMean(internal_computeMfccs(f_audio, samplingFreq,parameter),2);

numTotalFilters = parameter.MFCC_numLinearFilters+parameter.MFCC_numLogFilters;

% Now we find the desired magnitude response for the filter we are going to design
% by minimizing a least squares distance between the mfcc_dest and the mean MFCC
% computed from the filtered audio. This is a two step process.

% Step1: First we have to reconstruct the required gain/attenuation on mel
% spectrum from the information we have.
% That means we blow up a 12-13 dimensional vector back to 40, and
% hence there are many possible solutions. One way to find a solution is to
% expand the coeefficient vector to 40 dimensions by setting all unknown
% entries to 0, and then apply an inverse DCT, which is the transposed
% matrix (orthogonal). Another solution is to find the smoothest solution
% to the problem DCTcropped * x == f_mfccDestMean-f_mfccMean, and we get
% that by solving the corresponding tikhonov regularised least squares
% problem: ||DCTcropped * x - f_mfccDestMean-f_mfccMean||^2_2 +
% sigma*||x||_2 for some sigma>0. However, as it turns out, the direct
% inversion of the DCT is actually yielding more or less exactly this
% solution. Therefore, we don't have to do much about it here usually. Only
% if one prefers an even smoother solution (at the cost of not htting
% DCTcropped * x == f_mfccDestMean-f_mfccMean exactely), then the tikhonov
% solution might be helpful, so I'll leave the code in for now.
useTikhonovInversionStep1 = 0;
if useTikhonovInversionStep1
    DCTcropped = internal_getCroppedDCTMatrix(numTotalFilters,parameter.MFCC_cepstralCoeffsToKeep);
    
    %find the smoothest solution to this problem using a binary search on
    %candidates for the tikhonov lambda parameter.
    b = [(f_mfccDestMean-f_mfccMean);zeros(numTotalFilters,1)];
    lambda_candidates = sort(10.^-[0:0.2:5]);
    idxLow = 1;
    idxHigh = length(lambda_candidates);
    while idxHigh-idxLow>1
        fprintf('Low: %d   High: %d\n',idxLow,idxHigh)
        lambda = lambda_candidates(round( (idxHigh+idxLow)/2 ));
        A = [DCTcropped;lambda*eye(numTotalFilters)];
        logL_smooth = A \ b;
        euclDistance = norm(f_mfccMean + DCTcropped*logL_smooth - f_mfccDestMean);
        if euclDistance > 10^-3
            idxHigh = round( (idxHigh+idxLow)/2 )-1;
        else
            idxLow = round( (idxHigh+idxLow)/2 );
        end
    end
    lambda = lambda_candidates(idxLow);
    disp(lambda);
else
    b = zeros(numTotalFilters,1);
    b(parameter.MFCC_cepstralCoeffsToKeep) = f_mfccDestMean-f_mfccMean;
    logL_smooth = dct(eye(numTotalFilters))' * b;
end

% step 2: Here, just as above, we could transpose the mel filterbank and
% use its analysis functions as synthesis function, and indeed this kind
% of works OK. However, the mel filterbank matrix is not orthogonal, and
% when we would compute MFCCs from the reconstructed spectrum, we would get
% different MFCC. Hence we can't get around solving an actual system of
% equations here. Since this is underdetermined, we apply tikhonov
% regularization to get the smoothest solution. We also have to account for
% the preamphasis filter that was potentially applied
melFilterWeights = internal_getMelFilterbank(sideinfo.fftSize,samplingFreq,parameter);
% melFilterWeightsOrg = melFilterWeights;
if parameter.MFCC_applyPreemphasis
    freqRespPreemphasisFilter = abs(fft(parameter.MFCC_preemphasisFilter,sideinfo.fftSize));
    freqRespPreemphasisFilter = freqRespPreemphasisFilter(1:sideinfo.fftSize/2+1);
    melFilterWeights = melFilterWeights * diag(freqRespPreemphasisFilter);
end

numFreqBins = sideinfo.fftSize/2+1;
numTikhonovEntries = sideinfo.fftSize/2+1;
mainA = sparse([],[],[],numTotalFilters*numFreqBins,numFreqBins);
for f=1:numFreqBins
    mainA((f-1)*numTotalFilters+1:f*numTotalFilters,f) = melFilterWeights(:,f);
end
b = (diag(10.^logL_smooth) * melFilterWeights);
b = [b(:);zeros(numTikhonovEntries,1)]; % all cols are stacked to form one long upright vector
lambda_candidates = parameter.TikhonovLambdaCandidates;
idxLow = 1;
idxHigh = length(lambda_candidates);
runOnce = 0;
fullL_smooth = [];
while (~runOnce) || (idxHigh-idxLow>1)
    runOnce = 1;
    candIdx = round( (idxHigh+idxLow)/2 );
    lambda = lambda_candidates(candIdx);
    A = [mainA;lambda*eye(numTikhonovEntries)];
    fullL_smoothTemp = A \ b;
    if any(fullL_smoothTemp<0)
        %Usually, this does not happen. This is just in case. lsqnonneg is much slower
        warning('lsqnonneg is used to solve the Tikhonov reg. LS problem. This is just a comment to let you know that the method is slower than usual..')
        fullL_smoothTemp = lsqnonneg(A,b);
    end
    euclDistance = sqrt(sum(sum(  (diag(10.^logL_smooth)*melFilterWeights - melFilterWeights*diag(fullL_smoothTemp)).^2)));
    if euclDistance > parameter.TikhonovMaxError
        idxHigh = candIdx-1;
    else
        idxLow = candIdx;
        fullL_smooth = fullL_smoothTemp;
    end
end
if isempty(fullL_smooth)
    lambda = min(lambda_candidates);
    A = [mainA;lambda*eye(numTikhonovEntries)];
    fullL_smooth = lsqnonneg(A,b);
end

% now we 'scale' (subtract in log scale) the frequency response globally to
% get the least amount of overall gain. Here, the idea is that global
% scaling is just useful to get the same first MFCC coefficient as in the
% destination recording. However, since we normalize the recoring in the
% end anyway, that would be lost anyway. Since we don't have information
% about what to do with areas that are not covered by MFCC analysis
% function, we just try to keep the gain/attenuation everywhere as low as
% possible to get the remaining MFCCs done and don't do anything in the
% areas without support from the analysis functions. That gives room to
% design a more simple filter afterwards.
fullL_smooth_dB = zeros(size(fullL_smooth));
maxAtt_db = -20;
meanScaler = mean(20*log10(fullL_smooth(abs(fullL_smooth)>10^(maxAtt_db/20))));
fullL_smooth_dB(abs(fullL_smooth)>10^(maxAtt_db/20)) = 20*log10(fullL_smooth(abs(fullL_smooth)>10^(maxAtt_db/20))) - meanScaler;

% design filter (fir2 is linear phase)
filter_magFreqResp_dB_dest = fullL_smooth_dB;
filter_magFreqResp_linear_dest = 10 .^ (filter_magFreqResp_dB_dest/20);
h = fir2(parameter.filterOrder,sideinfo.freqinfoFFT/(samplingFreq/2),filter_magFreqResp_linear_dest);
   
% apply filter
parameterApplyImpulseResponse.loadInternalIR = 0;
parameterApplyImpulseResponse.impulseResponse = h;
parameterApplyImpulseResponse.impulseResponseSampFreq = samplingFreq;
parameterApplyImpulseResponse.normalizeOutputAudio = 1;
parameterApplyImpulseResponse.averageGroupDelayOfFilter = round(parameter.filterOrder/2);
[f_audio_out,timepositions_afterDegr] = degradationUnit_applyImpulseResponse(f_audio, samplingFreq, timepositions_beforeDegr, parameterApplyImpulseResponse);

if parameter.visualizations
    filterFreqResp_db = 20*log10(abs(fft(h,sideinfo.fftSize)));
    figure;
    plot([0:sideinfo.fftSize/2]/sideinfo.fftSize*samplingFreq,fullL_smooth_dB,'k');
    hold on
    plot([0:sideinfo.fftSize/2]/sideinfo.fftSize*samplingFreq,filterFreqResp_db(1:sideinfo.fftSize/2+1),'r');
    legend('required magnitude response','actual magnitude response of filter')
    
    parameter.MFCC_additionalFilter = 10.^(fullL_smooth_dB/20);
    f_mfccOutputMeanTemp = internal_infNanRobustMatrixMean(internal_computeMfccs(f_audio, samplingFreq,parameter),2);
    parameter.MFCC_additionalFilter = [];
    f_mfccOutputMean = internal_infNanRobustMatrixMean(internal_computeMfccs(f_audio_out, samplingFreq,parameter),2);
    figure;
    plot(f_mfccMean,'b');
    hold on
    plot(f_mfccDestMean,'k');
    plot(f_mfccOutputMeanTemp,'g');
    plot(f_mfccOutputMean,'r');
    legend('average MFCC of input audio','average MFCC of destination', 'average MFCC of input audio after filtering in freq domain','average MFCC of input audio after filtering in time domain')
end

end


function [f_logMelSpectra,sideinfo] = internal_computeLogMelSpectra(f_audio,fs,parameter)
% the code here was essentially taken from Malcalm Slaney's auditory
% toolbox and then simplified/rewritten where possible

% setting some parameters
windowSize =  2^nextpow2(fs*parameter.MFCC_windowSizeSec);
fftSize = 2*windowSize;
featuresPerSecond = parameter.MFCC_featuresPerSecond;
applyPreemphasis = parameter.MFCC_applyPreemphasis;
preemphasisFilter = parameter.MFCC_preemphasisFilter;
additionalFilter = parameter.MFCC_additionalFilter;

[melFilterWeights,sideinfoMel] = internal_getMelFilterbank(fftSize,fs,parameter);

% Filter the input with the preemphasis filter.  Also figure how
% many columns of data we will end up with.
preEmphasized = mean(f_audio,2);
if applyPreemphasis
    preEmphasized = filter(preemphasisFilter, 1, preEmphasized);  %oh boy. This is too weak for high sampling rates...
end
nOverlap = windowSize - round(fs/featuresPerSecond);

% spectrogram
[f_spec,freqinfo,timeinfo] = spectrogram(preEmphasized,hamming(windowSize),nOverlap,fftSize,fs);
if ~isempty(additionalFilter)
    f_spec = f_spec.* repmat(additionalFilter(:),1,size(f_spec,2));
end
f_logMelSpectra = log10(melFilterWeights * abs(f_spec));

if nargout > 1
    sideinfo.timeinfo = timeinfo;
    sideinfo.freqinfoFFT = freqinfo;
    sideinfo.allFreqs = sideinfoMel.freqs;
    sideinfo.centreFrequencies = sideinfoMel.center;
    sideinfo.lowerFrequencies = sideinfoMel.lower;
    sideinfo.upperFrequencies = sideinfoMel.upper;
    sideinfo.numTotalFilters = parameter.MFCC_numLinearFilters + parameter.MFCC_numLogFilters;
    sideinfo.windowSize = windowSize;
    sideinfo.fftSize = fftSize;
end

end

function DCTcropped = internal_getCroppedDCTMatrix(sizeDCT,cepstralCoeffsToKeep)
DCTcropped = dct(eye(sizeDCT));
DCTcropped = DCTcropped(cepstralCoeffsToKeep,:);
end

function [melFilterWeights,sideinfo] = internal_getMelFilterbank(fftSize,fs,parameter)
lowestFrequency = parameter.MFCC_lowestFrequencyHz;
numLinearFilters = parameter.MFCC_numLinearFilters;
numLogFilters = parameter.MFCC_numLogFilters;
numTotalFilters = numLinearFilters + numLogFilters;
linearFreqSpacingHz = parameter.MFCC_linearFreqSpacingHz;
logFreqSpacing = parameter.MFCC_logFreqSpacing;

% Now figure the band edges.  Interesting frequencies are spaced
% by linearSpacing for a while, then go logarithmic.  First figure
% all the interesting frequencies.  Lower, center, and upper band
% edges are all consequtive interesting frequencies.
freqs = lowestFrequency + (0:numLinearFilters-1)*linearFreqSpacingHz;
freqs(numLinearFilters+1:numTotalFilters+2) = ...
    freqs(numLinearFilters) * logFreqSpacing.^(1:numLogFilters+2);

lower = freqs(1:numTotalFilters);
center = freqs(2:numTotalFilters+1);
upper = freqs(3:numTotalFilters+2);

% We now want to combine FFT bins so that each filter has unit
% weight, assuming a triangular weighting function.  First figure
% out the height of the triangle, then we can figure out each
% frequencies contribution
melFilterWeights = zeros(numTotalFilters,fftSize/2+1);
triangleHeight = 2./(upper-lower);
fftFreqs = (0:fftSize/2)/fftSize*fs;

for chan=1:numTotalFilters
    melFilterWeights(chan,:) = ...
        (fftFreqs > lower(chan) & fftFreqs <= center(chan)).* ...
        triangleHeight(chan).*(fftFreqs-lower(chan))/(center(chan)-lower(chan)) + ...
        (fftFreqs > center(chan) & fftFreqs < upper(chan)).* ...
        triangleHeight(chan).*(upper(chan)-fftFreqs)/(upper(chan)-center(chan));
end

sideinfo.freqs = freqs;
sideinfo.lower = lower;
sideinfo.center = center;
sideinfo.upper = upper;

end

function [f_mfcc,sideinfo] = internal_computeMfccsFromLogMelSpectra(f_logMelSpectra,parameter)

cepstralCoeffsToKeep = parameter.MFCC_cepstralCoeffsToKeep;

% The following is equivalent to the auditory toolbox code, just shorter
% and more versatile.
DCTcropped = internal_getCroppedDCTMatrix(size(f_logMelSpectra,1),cepstralCoeffsToKeep);

f_mfcc = DCTcropped * f_logMelSpectra;

if nargout>1
    sideinfo.DCTcropped = DCTcropped;
end

end

function [f_mfcc,sideinfo] = internal_computeMfccs(f_audio,fs,parameter)

if nargout>1
    [f_logMelSpectra,sideinfo] = internal_computeLogMelSpectra(f_audio,fs,parameter);
    [f_mfcc,sideinfo2] = internal_computeMfccsFromLogMelSpectra(f_logMelSpectra,parameter);
    sideinfo.f_logMelSpectra = f_logMelSpectra;
    sideinfo.DCTcropped = sideinfo2.DCTcropped;
else
    f_logMelSpectra = internal_computeLogMelSpectra(f_audio,fs,parameter);
    f_mfcc = internal_computeMfccsFromLogMelSpectra(f_logMelSpectra,parameter);
end

end

function [output,locationWoutProblems] = internal_infNanRobustMatrixMean(matrix,dim)

summedMatrix = sum(matrix,mod(dim,2)+1); % 'preserves' nan and inf
locationWoutProblems = ~isinf(summedMatrix) & ~isnan(summedMatrix);

if dim == 1
    output = mean(matrix(locationWoutProblems,:),dim);
elseif dim == 2
    output = mean(matrix(:,locationWoutProblems),dim);
end

end







