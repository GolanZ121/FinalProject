clc; clear; close all;

%% Xor Mask & Cyclic Buffer
Fs = 15.36e6;
start = 4149;
buffer_size = 4236;

filename = 'raw_bits-20260423T161106Z-3-001\raw_bits/raw_bits_12_Feb_2026_09_30_39_442.txt';
fileID = fopen(filename,'r');
formatSpec = '%s';
bits = fscanf(fileID,formatSpec);
fclose(fileID);
bits = dec2bin(hex2dec(reshape(bits, numel(bits), [])))';
bits = bits(:);

mask = gen_gold_code;
demasked_bits = xor(bits,mask);

bits_after_cycle = cyclic_buffer(demasked_bits,start,buffer_size);

%% ECC & interliver (turbo decoder)

ProccesedBits=bits_after_cycle;
interliver= mod((43 * (1/1408) + 88 * (1/1408).^2), 1408);
ConvCodeRate=1/2;
SubBlockLength=1412;
DummyFormat=[1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31,2,18,10,26,6,22,14,30,4,20,12,28,8,24,16,32];
DummySize=28;
DummyMatSize=[45, 32];

data=TurboDecoder(ProccesedBits, SubBlockLength, interliver, DummySize, DummyMatSize, DummyFormat);

%% CRC

deECCbits = data;

deECCbitsLen=1408;
CRCLen=24; %at the end of the data
polynom=[1 1 0 1 1 1 1 1 0 0 1 1 0 0 1 0 0 1 1 0 0 0 0 1 1]; % 0 --> n=24 
init = zeros(1,CRCLen);

IsValidCRC = CRC_Vailidation(deECCbits, polynom, init);

%% Xor Mask & Cyclic Buffer functions

% this function creates a xor mask 
function seq = gen_gold_code(Nc, L, seed)
arguments
    Nc (1,1) = 1600
    L (1,1) = 7200
    seed (1,1) = 0x12345678
end

x1 = zeros(Nc + L + 31, 1);
x2 = zeros(Nc + L + 31, 1);
x2(1:32) = flip(dec2bin(seed, 32));
x1(1) = 1;
for n = 1:(Nc + L)
    x1(n + 31) = xor(x1(n + 3), x1(n));
    x2(n + 31) = xor(xor(x2(n + 3), x2(n + 2)), xor(x2(n+1), x2(n)));
end
seq = xor(x1(Nc + 1:Nc+L), x2(Nc + 1:Nc+L));
end

% this function removes a cyclic buffer
function bits_after_cycle = cyclic_buffer(bits,start,buffer_size)
% bits - bits of info
% start - start of info
% size - size of info
stop = buffer_size - (length(bits) - start); % where in the bits array to stop after the cyclic copy
bits_after_cycle = [bits(start:end);bits(1:stop - 1)];
end

%% Turbo Decoder functions

function data=TurboDecoder(ProccesedBits, SubBlockLength, interliver, DummySize, DummyMatSize, DummyFormat)

% input: parameters and processed (after xor and cyclic shift removal) bits
% output: data bits for CRC 

[S1, P1, P2]= PreDecoderDecode(ProccesedBits, SubBlockLength, interliver, DummySize, DummyMatSize, DummyFormat);

ConstraintLen= 4; % because the polynomial degree is 3 
FeedBackConn =13; % because 13 in octal base is 11 which is 1011 in binary base which is equal to 1+D^2+D^3 which is the polynom of the encryption
CodeGen= [15 13]; % 15[oct]= 13[dec]= 1101[bin]=1+D+D^3 (G)
trellis = poly2trellis(ConstraintLen,CodeGen,FeedBackConn); 

tbdepthPres=[1/2 5;2/3 7.5;3/4 10;5/6 15]; % common tbdepth by CodeRate according to matlab documentation
tbdepth=tbdepthPres(tbdepthPres(:,1)==ConvCodeRate,2)*(ConstraintLen-1); %because rate 1/2 code has a traceback depth of 5 × (ConstraintLength – 1)
opmode='trunc'; % because the encoder is assumed to have started at the all-zeros state.
dectype='hard'; % because the input values are 0 or 1.

S2= vitdec([S1,P1],trellis,tbdepth,opmode,dectype);
data= vitdec([S2,P2],trellis,tbdepth,opmode,dectype);

end

function [S, P1, P2]=PreDecoderDecode(ProccesedBits, SubBlockLength, interliver, DummySize, DummyMatSize, DummyFormat)

InterleavedS=ProccesedBits(1:SubBlockLength);
InterleavedInerlacedP=ProccesedBits(SubBlockLength+1:end);

[InterleavedP1, InterleavedP2] =Deinterlacer(InterleavedInerlacedP);
DeDummiedS=DeDummy(InterleavedS, DummySize, DummyMatSize, DummyFormat);
DeDummiedP1=DeDummy(InterleavedP1, DummySize, DummyMatSize, DummyFormat);
DeDummiedP2=DeDummy(InterleavedP2, DummySize, DummyMatSize, DummyFormat);

S=DeInterliver(interliver, DeDummiedS);
P1=DeInterliver(interliver, DeDummiedP1);
P2=DeInterliver(interliver, DeDummiedP2);

end

function DeInterlivedSubBlock=DeInterliver(interliver, InterleavedSubBlock)

DeInterlivedSubBlock = deintrlv(InterleavedSubBlock, interliver);

DeInterlivedSubBlock=DeInterlivedSubBlock(1:end-4);

end

function DeDummiedBlock=DeDummy(InterleavedSubBlock, DummySize, DummyMatSize, DummyFormat)

%reshape to mat by columns 
ColBlock=reshape(InterleavedSubBlock,DummyMatSize);

%insert dummies
DummiBlock=[NaN(DummyMatSize(1),1),ColBlock];

%reshape to row by rows
DummiBlock=reshape(DummiBlock.',1,[]);

%remove unneeded dummies
DummyJumps=(DummySize-DummyMatSize(2))*2; % (32-28)*2=8, might be a coincidence
DummiBlock(1:DummyMatSize(1)*DummyJumps:end)=[];

%reshape to mat by rows
DummiBlock=reshape(DummiBlock,flip(DummyMatSize)).';

%mix back the columns
DummiBlock=DummiBlock(:,DummyFormat);

%reshape to row by columns
DummiBlock=reshape(DummiBlock,1,[]);

%remove dummies
DeDummiedBlock=DummiBlock(DummySize+1:end);

end

function [InterleavedP1, InterleavedP2] = Deinterlacer(InterleavedInerlacedP)
InterleavedP1= InterleavedInerlacedP(1:2:end);
InterleavedP2= InterleavedInerlacedP(2:2:end);
end

%% CRC functions

function IsValidCRC = CRC_Vailidation(deECCbits, polynom, init)

% validates the check value of the Cyclic redundancy check 
% parameters:
% deECCbits - the bits after decoding the xor mask, cyclic buffer and error correction code (interliver)
% polynom - the characteristic polynomial 
% init - the initialization (sequence that was concated to the data)

[~,res] = polydiv(deECCbits,polynom);

IsValidCRC = (res == init);

end
