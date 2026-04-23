%% Scripts to examine the samples
clc; clear; close all;
% Params
Fs = 10e6;

% load the raw_samples file 
filename = 'raw_samples_12_Feb_2026_09_30_40_676_fs_10MHz.32fc';
fid = fopen(filename, 'r');

% Read all data as 32-bit floats
data = fread(fid, inf, 'float32'); 
fclose(fid);
numSamples = length(data);
TotalTime = numSamples / Fs;
% Reshape data into a time series format

% time = (0:numSamples-1) / Fs;
% plot(time, [real(data), imag(data)]);
% legend("Real", "Imaginary");

% nOverlap = 0
s = spectrogram(data)
% [s,f,t] = spectrogram(data,win,nOverlap,freqSpec,Fs)