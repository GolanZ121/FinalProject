clc; clear; close all;

%% definisions
Fs=10e6;%Hz (for sdr60)
BB=2.4e9;%Hz
RecordTime=3.2;%seconds - for sdr
PacketTime=643e-6;%seconds

%% Scripts to examine the samples

% load the raw_samples file 
filename = 'raw_samples/raw_samples_12_Feb_2026_09_30_39_442_fs_10MHz.32fc';
fid = fopen(filename, 'r');

% Read all data as 32-bit floats
data = fread(fid, [2, Inf], 'float32=>double');
fclose(fid);
data = complex(data(1,:), data(2,:));

numSamples = length(data);
TotalTime = numSamples / Fs;
SampleTime=1/Fs;

time = (0:numSamples-1) / Fs;
Symbol4Time =(59.62/SampleTime:59.69/SampleTime)/ Fs;
Symbol6Time = (59.76/SampleTime:59.83/SampleTime)/ Fs;
Symbol4data=data(59.62/SampleTime:59.69/SampleTime);
Symbol6data = data(59.76/SampleTime:59.83/SampleTime);

figure;
plot(real(Symbol4data));
hold on;
plot(imag(Symbol4data));




































































