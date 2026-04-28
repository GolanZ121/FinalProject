clc; clear; close all;

%% definitions

% deECCbits = []; % MUST IMPLEMENT !!!
deECCbits = zeros(1,1408);

% CRC_CheckValue parameters
deECCbitsLen=1408;
CRCLen=24; %at the end of the data
polynom=[1 1 0 1 1 1 1 1 0 0 1 1 0 0 1 0 0 1 1 0 0 0 0 1 1]; % 0 --> n=24 
init = zeros(CRCLen,1).';

%% test
IsValidCRC = CRC_Vailidation(deECCbits, polynom, init);

%% CRC functions
function CheckValue = CRC_CheckValue(deECCbits, CRCLen, deECCbitsLen, polynom,init)

% finds the check value of the Cyclic redundancy check 
% parameters:
% deECCbits - the bits after decoding the xor mask, cyclic buffer and error correction code (interliver)
% CRCLen - the length of the CRC sequence (n)
% deECCbitsLen - the length of the bit sequence (data+CRC_CheckValue)
% polynom - the characteristic polynomial 

bits=[deECCbits(1:deECCbitsLen-CRCLen),init];

%
% for i=1:deECCbitsLen
%     deECCbits_xor_polynom=xor(bits(i:i+CRCLen), polynom);
%     bits=[deECCbits_xor_polynom,bits(i+CRCLen+1:end)];
% end
[bits,~] = polydiv(bits,polynom);

CheckValue=bits(deECCbitsLen-CRCLen:end);

end

function IsValidCRC = CRC_Vailidation2(deECCbits, CRCLen, deECCbitsLen, polynom,init)

% validates the check value of the Cyclic redundancy check 
% parameters:
% deECCbits - the bits after decoding the xor mask, cyclic buffer and error correction code (interliver)
% CRCLen - the length of the CRC sequence (n)
% deECCbitsLen - the length of the bit sequence (data+CRC_CheckValue)
% polynom - the characteristic polynomial 

CheckValue = CRC_CheckValue(deECCbits, CRCLen, deECCbitsLen, polynom,init);

bits=[deECCbits(1:deECCbitsLen-CRCLen),CheckValue];

% for i=1:deECCbitsLen
%     deECCbits_xor_polynom=xor(bits(i:i+CRCLen), polynom);
%     bits=[deECCbits_xor_polynom,bits(i+CRCLen+1:end)];
% end
[bits,res] = polydiv(bits,polynom);

IsValidCRC = (bits(deECCbitsLen-CRCLen:end)== zeros(CRCLen,1).');

end

function IsValidCRC = CRC_Vailidation(deECCbits, polynom, init)

% validates the check value of the Cyclic redundancy check 
% parameters:
% deECCbits - the bits after decoding the xor mask, cyclic buffer and error correction code (interliver)
% polynom - the characteristic polynomial 
% init - the initialization (sequence that was concated to the data)

[~,res] = polydiv(deECCbits,polynom);

IsValidCRC = (res == init);

end


function success = CRC_test(deECCbits, CRCLen, deECCbitsLen, polynom, init, CheckValue)

% validates the check value of the Cyclic redundancy check 
% parameters:
% deECCbits - the bits after decoding the xor mask, cyclic buffer and error correction code (interliver)
% CRCLen - the length of the CRC sequence (n)
% deECCbitsLen - the length of the bit sequence (data+CRC_CheckValue)
% polynom - the characteristic polynomial 

bits=[deECCbits(1:deECCbitsLen-CRCLen),init];

% for i=1:deECCbitsLen
%     deECCbits_xor_polynom=xor(bits(i:i+CRCLen), polynom);
%     bits=[deECCbits_xor_polynom,bits(i+CRCLen+1:end)];
% end

[deECCbits,~] = polydiv(bits,polynom);

success = CRC_Vailidation(deECCbits, CRCLen, deECCbitsLen, polynom,init);

% success=(bits(deECCbitsLen-CRCLen:end)==CheckValue);

end

