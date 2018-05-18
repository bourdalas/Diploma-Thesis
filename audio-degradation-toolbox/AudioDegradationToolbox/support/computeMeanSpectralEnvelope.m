function computeMeanSpectralEnvelope(audiofile, fftLength)
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

[f_audio,samplingFreq]=wavread(audiofile);

% compute mean spectral vector for f_audio
[destMagFreqResp, destMagFreqResp_freqs] = internal_computeMeanSpectralVector(f_audio,samplingFreq,fftLength);
destMagFreqResp = internal_standardizeMagFreqResp(destMagFreqResp);

outputFilename = [audiofile(1:end-4),'_specEnv'];
save(outputFilename,'destMagFreqResp','destMagFreqResp_freqs')

end


function [f_meanmagspec_db,freqs] = internal_computeMeanSpectralVector(f_audio,fs,fftLength)

f_audio = mean(f_audio,2);

[f_spec,freqs,time] = spectrogram(f_audio,hanning(fftLength),fftLength/2,fftLength,fs);

f_magspec_db = 20 * log10(abs(f_spec));

f_magspec_db(:,isinf(sum(abs(f_magspec_db),1))) = []; % ignore rows with -inf/inf entries

f_meanmagspec_db = mean(f_magspec_db,2);

end

function magFreqResp = internal_standardizeMagFreqResp(magFreqResp)

maxRobust = max(magFreqResp(~isinf(magFreqResp)));

magFreqResp = magFreqResp - maxRobust;

magFreqResp(magFreqResp > 0) = 0;  % remaining positive inf
magFreqResp(magFreqResp < -80) = -80;  % remaining positive inf

end

