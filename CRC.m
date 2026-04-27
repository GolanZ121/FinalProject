% בעיה: אי אפשר לבדוק את החלק הזה בלי לקבל פקטה חתוכה ולדעת מה אנחנו
% מצפים לקבל וגם אני לא מצליחה לפתוח את הקובץ כביטים (אני פשוט לא בטוחה שפתחתי אותו תקין)
clc; clear; close all;

%% definitions
Fs = 10e6;

% load the raw_samples file 
filename = 'C:\Users\ZBOOK\FinalProject\processed_bits\processed_bits12_Feb_2026_09_30_39_442.txt';
fid = fopen(filename, 'r');

% Read all data
deECCbits = dec2bin(fread(fid));
fclose(fid);

% CRC_CheckValue parameters
deECCbitsaLen=1408;
CRCLen=24; %at the end of the data
polynom=[1 1 0 1 1 1 1 1 0 0 1 1 0 0 1 0 0 1 1 0 0 0 0 1 1]; % 0 --> n=24 
init=0;

%% test
success = CRC_test(deECCbits, CRCLen, deECCbitsaLen, polynom, zeros(CRCLen,1).');


%% CRC
function CheckValue = CRC_CheckValue(deECCbits, CRCLen, deECCbitsaLen, polynom)

% finds the check value of the Cyclic redundancy check 
% parameters:
% deECCbits - the bits after decoding the xor mask, cyclic buffer and error correction code (interliver)
% CRCLen - the length of the CRC sequence (n)
% deECCbitsaLen - the length of the bit sequence (data+CRC_CheckValue)
% polynom - the characteristic polynomial 

bits=[deECCbits(1:deECCbitsaLen-CRCLen),zeros(CRCLen,1).'];

%
for i=1:deECCbitsaLen
    deECCbits_xor_polynom=xor(bits(i:i+CRCLen), polynom);
    bits=[deECCbits_xor_polynom,bits(i+CRCLen+1:end)];
end

CheckValue=bits(deECCbitsaLen-CRCLen:end);

end

function IsValidCRC = CRC_Vailidation(deECCbits, CRCLen, deECCbitsaLen, polynom)

% validates the check value of the Cyclic redundancy check 
% parameters:
% deECCbits - the bits after decoding the xor mask, cyclic buffer and error correction code (interliver)
% CRCLen - the length of the CRC sequence (n)
% deECCbitsaLen - the length of the bit sequence (data+CRC_CheckValue)
% polynom - the characteristic polynomial 

CheckValue = CRC_CheckValue(deECCbits, CRCLen, deECCbitsaLen, polynom);

bits=[deECCbits(1:deECCbitsaLen-CRCLen),CheckValue];

for i=1:deECCbitsaLen
    deECCbits_xor_polynom=xor(bits(i:i+CRCLen), polynom);
    bits=[deECCbits_xor_polynom,bits(i+CRCLen+1:end)];
end

IsValidCRC=bits(deECCbitsaLen-CRCLen:end)==zeros(CRCLen,1).';

end

function success = CRC_test(deECCbits, CRCLen, deECCbitsaLen, polynom, CheckValue)

% validates the check value of the Cyclic redundancy check 
% parameters:
% deECCbits - the bits after decoding the xor mask, cyclic buffer and error correction code (interliver)
% CRCLen - the length of the CRC sequence (n)
% deECCbitsaLen - the length of the bit sequence (data+CRC_CheckValue)
% polynom - the characteristic polynomial 

bits=[deECCbits(1:deECCbitsaLen-CRCLen),zeros(CRCLen,1).'];

for i=1:deECCbitsaLen
    deECCbits_xor_polynom=xor(bits(i:i+CRCLen), polynom);
    bits=[deECCbits_xor_polynom,bits(i+CRCLen+1:end)];
end

success = CRC_Vailidation(deECCbits, CRCLen, deECCbitsaLen, polynom);

% success=(bits(deECCbitsaLen-CRCLen:end)==CheckValue);

end