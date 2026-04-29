%% TEST TO DELETE LATER !!!
filename = 'C:\Users\ZBOOK\FinalProject\processed_bits\processed_bits12_Feb_2026_09_30_39_442.txt';
fileID = fopen(filename,'r');
formatSpec = '%s';
bits = fscanf(fileID,formatSpec);
fclose(fileID);
bits = dec2bin(hex2dec(reshape(bits, numel(bits), [])))';
bits = bits(:).';
