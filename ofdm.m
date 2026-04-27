%% clean
clc; clear; close all;
%% Scripts to examine the samples
% Params
Fs = 15.36e6;

% load the raw_samples file 
filename = 'synced_samples-20260423T161109Z-3-001/synced_samples/synced_samples_12_Feb_2026_09_30_39_442_fs_15.36MHz.32fc';
fileID = fopen(filename, 'r');
data = fread(fileID, [2, Inf], 'single');
data = data(1,:) + 1i*data(2,:);
fclose(fileID);

numSamples = length(data);
TotalTime = numSamples / Fs;
% Reshape data into a time series format
time = (0:numSamples-1) / Fs;
plot(time, data);
legend("Real", "Imaginary");

%% Script to examine ofdm function
% Params
Ncp = 72;
ex_Ncp = 80;
Nsc = 1024;
start_sc = 213;
stop_sc = 813;
num_of_symbols = 9;
NFFT = 1024;

%% FIX ONE BLOCK: 

% block2 = data(81:(1024+80));
% x = fftshift(fft(block2));


start_block4 = 1024 * 4 + 80 + 72 * 4 + 1;
start_block6 = 1024 * 5 + 80 +72 * 4 + 1;
block4 = data(start_block4:(1024 + start_block4 - 1));
x = fftshift(fft(block4));
phases = linspace(0,2,360);

scatterplot(x)

%% CONTINUE
demodulated_data = ofdm_demod(data,Ncp,ex_Ncp,Nsc,start_sc,stop_sc, num_of_symbols);

% load the raw_samples file 
filename = 'raw_bits-20260423T161106Z-3-001\raw_bits/raw_bits_12_Feb_2026_09_30_39_442.txt';
fileID = fopen(filename,'r');
formatSpec = '%s';
bits = fscanf(fileID,formatSpec);
fclose(fileID);
binresult = dec2bin(hex2dec(reshape(bits, numel(bits)/2, [])));
%% functions
function data = ofdm_demod(packet,Ncp,ex_Ncp,Nsc,start_sc,stop_sc, num_of_symbols)
    % packet - the packet of information
    % Ncp - cyclic prefix, a partial copy of the end of the symbol to create  a
    % correlation peak at each start of symbol
    % ex_Ncp - extended cp
    % start_sc - is the start of the active subcarrier range
    % stop_sc - is the start of the active subcarrier range
    % num_of_symbols - a given information about num of symbols in a packet
    
    clean_packet = remove_cp(packet,Ncp,ex_Ncp,Nsc,num_of_symbols); % packet without cp
    
    slot = reshape(clean_packet,Nsc,[]); % create the slot matrix
    slot_fft = fftshift(fft(slot))./(sqrt(Nsc)); % according to the formula (not sure if needed in our case)
    slot_fft = slot_fft(start_sc:stop_sc,:); % cut the irrelevant frequencies (frequency domain)
    slot_fft = [slot_fft(1:600,:),slot_fft(602:end,:)]; % the mid subcarrier is null
    slot_fft = [slot_fft(:,2:3),slot_fft(:,5),slot_fft(:,7:end)];

    % transform the slot to bits of data
    data = slot_fft(:);
    data = pskdemod(data,4);
    data = dec2bin(data)';
    data = data(:);
end

% this function receives a packet with cp and returns it without
function clean_packet = remove_cp(packet,Ncp,ex_Ncp,Nsc,num_of_symbols)
    % packet - the packet of information
    % Ncp - cyclic prefix, a partial copy of the end of the symbol to create  a
    % correlation peak at each start of symbol
    % ex_Ncp - extended cp
    % num_of_symbols - a given information about num of symbols in a packet

    cp_diffs = ex_Ncp - Ncp;  % We will use it to remove the extended cp later
    clean_packet = [];
    symbol_index = Ncp + 1;
    % at each iteration of the loop I cut the cp and paste the data
    for i = [1:num_of_symbols]
        %symbols 1 and 9  have extended cp
        if i == 1 || i == num_of_symbols
            symbol_index = symbol_index + cp_diffs ;
        end
        end_symbol_index = symbol_index + Nsc  - 1;
       clean_packet = [clean_packet,packet(symbol_index :end_symbol_index )] ;
       symbol_index = symbol_index + Nsc + Ncp; % calculate the next symbol index
    end
end