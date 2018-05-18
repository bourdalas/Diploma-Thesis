# Diploma-Thesis
Music Information Retrieval System Implementation in Matlab for Unix.

The project includes:

Audio Fingerprinting Toolbox for Music Identification Retrival System: It uses two encoding techinques for wav files, that use as basic tools the Short Time Fourier Transform, and Constant-Q Transformation accordingly.

Audio Degradation Toolbox: Can be used in order to distort the input signals, so that the efficincy of the fingerprint is tested.
https://code.soundsoftware.ac.uk/projects/audio-degradation-toolbox

CQT Toolbox: Used in order to compute the Constant-Q Transformation
https://www.cs.tut.fi/sgn/arg/CQT/

Matlab-LMDB-Wrapper: Used in order to use LMDB as the system's database engine. 

https://github.com/kyamagu/matlab-lmdb

In order for the system to work, the user must build the Matlab Wrapper for LMDB database
Instructions:
  Specify the MATLABDIR path.
  make MATLABDIR=/usr/local/MATLAB/MATLAB_VERSION
  make test
