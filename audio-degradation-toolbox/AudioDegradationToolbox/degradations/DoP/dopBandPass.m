function degradation_config = dopBandPass()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preset: dopBandPass
% Programmer: Matthias Mauch
%
% Description:
% Parametrisation of the Degradation Unit of the same name
%         
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

degradation_config(1).methodname = 'degradationUnit_applyLowpassFilter';
degradation_config(1).parameter.stopFrequency = 6000;
degradation_config(1).parameter.passFrequency = 4000;

degradation_config(2).methodname = 'degradationUnit_applyHighpassFilter';
degradation_config(2).parameter.stopFrequency = 66.66;
degradation_config(2).parameter.passFrequency = 200;

% degradation_config(3).methodname = 'degradationUnit_applyDynamicRangeCompression';

degradation_config(3).methodname = 'adthelper_normalizeAudio';