clc; clear; close all;

%% definitions
Fs = 10e6;

% load the raw_samples file 
filename = 'C:\Users\ZBOOK\FinalProject\raw_bits\raw_bits_12_Feb_2026_09_30_39_442.txt';
fid = fopen(filename, 'r');

% Read all data as 32-bit floats
deECCbits = fread(fid, [2, Inf], 'float32=>double');
fclose(fid);

deECCbitsaLen=1408;
CRCLen=24; %at the end of the data

%% CRC
function ProcessedBits = crc(deECCbits, CRCLen, deECCbitsaLen)

% Cyclic redundancy check decoding and removal
% parameters:
% deECbits - the bits after decoding the xor mask, cyclic buffer and error correction code (interliver)
% CRCLen - the length of the CRC sequence 
% deECCbitsaLen - the length of the bit sequence (data+CRC)



end