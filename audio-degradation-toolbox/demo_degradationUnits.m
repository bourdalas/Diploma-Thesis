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

% Note: Some degradations impose a delay/temporal distortion on the input
% data. The last example shows, given time positions before the
% degradation, how the corresponding time positions after the degradation
% can be retrieved

%%

clear
close all

addpath(genpath(fullfile(pwd,'AudioDegradationToolbox')));

pathOutputDemo = 'demoOutput/';
if ~exist(pathOutputDemo,'dir'), mkdir(pathOutputDemo); end

filenames = {
    'code/audio-degradation-toolbox/testdata/RWC_G39.wav';
    'code/audio-degradation-toolbox/testdata/RWC_G72.wav';
    'code/audio-degradation-toolbox/testdata/RWC_G84.wav';
    'code/audio-degradation-toolbox/testdata/RWC_P009m_drum.wav';
    'code/audio-degradation-toolbox/testdata/RWC-C08.wav';
    'code/audio-degradation-toolbox/testdata/session5-faure_elegie2c-001-0.wav';
    'code/audio-degradation-toolbox/testdata/175234__kenders2000__nonsense-sentence.wav';
    };

createSpectrograms = 1;

%%
% just copying original files to the demo folder
maxValueRangeVis = zeros(length(filenames));
for k=1:length(filenames)
    copyfile(filenames{k}, fullfile(pathOutputDemo,sprintf('00_Original_file%d.wav',k)))
    if createSpectrograms
        [f_audio,samplingFreq]=audioread(filenames{k});
        [s,f,t] = spectrogram(f_audio,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1)); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('00_Original_file%d.png',k)))
        maxValueRangeVis(k) = max(max(log10(abs(s)+1)));
    end
end

%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    % with default settings:
    %f_audio_out = degradationUnit_addNoise(f_audio, samplingFreq);
    
    % adjusting some parameters:
    parameter.snrRatio = 10; % in dB
    parameter.noiseColor = 'pink';  % convenient access to several noise types
    f_audio_out = degradationUnit_addNoise(f_audio, samplingFreq, [], parameter);
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_01_addNoise_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_01_addNoise_file%d.png',k)))    
    end
end

%%
% Some degradations delay the input signal. If some timepositions are given
% via timepositions_beforeDegr, the corresponding positions will be returned
% in timepositions_afterDegr.
timepositions_beforeDegr = [2, 3];

% Processing the time positions even works without processing the audio data (f_audio_out is empty afterwards)
parameter.loadInternalIR = 0;
parameter.impulseResponse = [1 -1 1 -1 1 -1 ] / 6;  % some impulse response
parameter.impulseResponseSampFreq = samplingFreq;
parameter.normalizeOutputAudio = 1;
[f_audio_out,timepositions_afterDegr] = degradationUnit_applyImpulseResponse([], [], timepositions_beforeDegr, parameter);

for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
        
    % time positions and audio can also be processed at the same time:
    [f_audio_out,timepositions_afterDegr] = degradationUnit_applyImpulseResponse(f_audio, samplingFreq, timepositions_beforeDegr, parameter);
    fprintf('degradation_applyFirFilter: adjusting time positions\n');
    for m=1:length(timepositions_afterDegr) fprintf('%g -> %g\n',timepositions_beforeDegr(m),timepositions_afterDegr(m)); end
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_02_applyImpulseResponse_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_02_applyImpulseResponse_file%d.png',k)))    
    end
end;

%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    parameter.snrRatio = 10; % in dB
    parameter.loadInternalSound = 1;
    parameter.internalSound = 'PubEnvironment1';
    f_audio_out = degradationUnit_addSound(f_audio, samplingFreq, [], parameter);
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_03_addSound_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_03_addSound_file%d.png',k)))    
    end
end;

%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    parameter.dsFrequency = 4000;
    f_audio_out = degradationUnit_applyAliasing(f_audio, samplingFreq, [], parameter);
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_04_applyAliasing_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_04_applyAliasing_file%d.png',k)))    
    end
end;

%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    parameter.percentOfSamples = 10;  % signal is scaled so that n% of samples clip
    f_audio_out = degradationUnit_applyClippingAlternative(f_audio, samplingFreq, [], parameter);
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_05_applyClipping_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_05_applyClipping_file%d.png',k)))    
    end
end;
%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    parameter.compressorSlope = 0.9;
    parameter.normalizeOutputAudio = 1;
    f_audio_out = degradationUnit_applyDynamicRangeCompression(f_audio, samplingFreq, [], parameter);
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_06_applyDynamicRangeCompression_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_06_applyDynamicRangeCompression_file%d.png',k)))    
    end
end;

%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    parameter.nApplications = 5;
    f_audio_out = degradationUnit_applyHarmonicDistortion(f_audio, samplingFreq, [], parameter);
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_07_applyHarmonicDistortion_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_07_applyHarmonicDistortion_file%d.png',k)))    
    end
end;

%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    parameter.LameOptions = '--preset cbr 32';
    f_audio_out = degradationUnit_applyMp3Compression(f_audio, samplingFreq, [], parameter);
    
    audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_08_applyMp3Compression_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_08_applyMp3Compression_file%d.png',k)))    
    end
end;

%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    parameter.changeInPercent = +5;
    f_audio_out = degradationUnit_applySpeedup(f_audio, samplingFreq, [], parameter);
    

%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_09_applySpeedup_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_09_applySpeedup_file%d.png',k)))    
    end
end;

%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    parameter.intensityOfChange = 3;
    parameter.frequencyOfChange = 0.5;
    f_audio_out = degradationUnit_applyWowResampling(f_audio, samplingFreq, [], parameter);
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_10_applyWowResampling_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_10_applyWowResampling_file%d.png',k)))    
    end
end;

%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    parameter.stopFrequency = 1000;
    f_audio_out = degradationUnit_applyHighpassFilter(f_audio, samplingFreq, [], parameter);
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_11_applyHighpassFilter_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_11_applyHighpassFilter_file%d.png',k)))    
    end
end;

%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    parameter.stopFrequency = 1000;
    f_audio_out = degradationUnit_applyLowpassFilter(f_audio, samplingFreq, [], parameter);
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_12_applyLowpassFilter_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_12_applyLowpassFilter_file%d.png',k)))    
    end
end;


%%
for k=1:length(filenames)
    [f_audio,samplingFreq]=audioread(filenames{k});
    
    % There are four ways to specify the destination mean spectrum: 1. by
    % loading example files provided with the toolbox, 2. by using specific
    % noise "color" profiles, 3. by providing the destination mean spectrum
    % using the parameter destMagFreqResp, 4. by providing audio data from
    % which the destination mean spectrum is computed.
    parameter.loadInternalMagFreqResp = 1; 
    parameter.internalMagFreqResp = 'Beethoven_Appasionata_Rwc';
    f_audio_out = degradationUnit_adaptiveEqualizer(f_audio, samplingFreq, [], parameter);
    
%     audiowrite(f_audio_out,samplingFreq,16,fullfile(pathOutputDemo,sprintf('Unit_13_adaptiveEqualizer_file%d.wav',k)));
    if createSpectrograms
        [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
        figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_13_adaptiveEqualizer_file%d.png',k)))
    end
end;

%%
if 1
    for k=1:length(filenames)
        [f_audio,samplingFreq]=audioread(filenames{k});
       
        [parameter.audioDataForDestMfcc,parameter.audioDataForDestMfcc_sf]=audioread('testdata/RWC_P009m_drum.wav');

        % After this unit, the mean MFCC vector of f_audio_out is almost identical to
        % the one of RWC_P009m_drum.wav
        parameter.visualizations = createSpectrograms;
        f_audio_out = degradationUnit_applyMfccMeanAdaption(f_audio, samplingFreq, [], parameter);
              
%         audiowrite('Unit_14_applyMfccMeanAdaption_file%d.wav',f_audio_out,samplingFreq,16,k)));
        if createSpectrograms
            [s,f,t] = spectrogram(f_audio_out,hamming(round(samplingFreq*0.093)),round(samplingFreq*0.093/2),[],samplingFreq);
            figure; imagesc(t,f,log10(abs(s)+1),[0 maxValueRangeVis(k)]); axis xy; colormap(hot); ylim([0,8000]); colorbar; print('-dpng', fullfile(pathOutputDemo,sprintf('Unit_14_applyMfccMeanAdaption_file%d.png',k)))
        end
    end
end



