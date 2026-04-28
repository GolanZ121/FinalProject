%% SDR side 
clc; clear; close all;
Fs = 10e6;
% Create USRP receiver
time_of_frame = 2*640e-6;
samples_per_frame = time_of_frame * Fs;

rx = comm.SDRuReceiver( ...
    'Platform', 'B210', ...
    'SerialNum', '356E63A', ...
    'CenterFrequency', 433e6, ...  
    'Gain', 30, ...
    'SampleRate', Fs, ...
    'SamplesPerFrame', samples_per_frame);

while true
    [data, len] = rx();
    spectrogram(double(data),100,80,100,Fs,'centered')
    title('Live B210 Signal');
    drawnow;
end