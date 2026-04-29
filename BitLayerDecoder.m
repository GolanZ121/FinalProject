
function BitLayerDecoder (bits)

% xor mask & cyclic buffer parameters
start = 4149;
buffer_size = 4236;

% xor mask
mask = GenGoldCode;
bits = de2bi(bin2dec(bits),'left-msb');
DemaskedBits = xor(bits,mask);

% cyclic mask
BitsAfterCyclicBufferRemoval = CyclicBufferRemover(DemaskedBits.',start,buffer_size);

% ECC & interlivers parameters (turbo decoder)

PreECCbits=BitsAfterCyclicBufferRemoval;
interliver= mod((43 * (1:1408) + 88 * (1:1408).^2), 1408);
ConvCodeRate=1/2;
SubBlockLength=1412;
DummyFormat=[1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31,2,18,10,26,6,22,14,30,4,20,12,28,8,24,16,32];
DummySize=28;
DummyMatSize=[45, 32];

% separate S, P1 and P2
InterleavedS=PreECCbits(1:SubBlockLength);
InterleavedInerlacedP=PreECCbits(SubBlockLength+1:end); 

[InterleavedP1, InterleavedP2] =Deinterlacer(InterleavedInerlacedP);

% dummy de-interliver
DummySpots=Dummy(SubBlockLength, DummySize, DummyMatSize, DummyFormat);

S=DeDummy(InterleavedS, DummySize, DummyMatSize, DummyFormat, DummySpots);
P1=DeDummy(InterleavedP1, DummySize, DummyMatSize, DummyFormat, DummySpots);
P2=DeDummy(InterleavedP2, DummySize, DummyMatSize, DummyFormat, DummySpots);

% viterbi and pi de-interliver
DecodedData=TurboDecoder(S, P1, P2, interliver, ConvCodeRate);


% CRC parameters
deECCbits = DecodedData;
deECCbitsLen=1408;
CRCLen=24; %at the end of the data
polynom=[1 1 0 1 1 1 1 1 0 0 1 1 0 0 1 0 0 1 1 0 0 0 0 1 1]; % 0 --> n=24 
init = zeros(1,CRCLen);

% CRC check
IsValidCRC = CRC_Vailidation(deECCbits, deECCbitsLen, polynom, init);

% parsing parameters
deECCbits=deECCbits(1:728);
deECCbits = dec2bin(deECCbits).';

% parsing
parsing(deECCbits)

end

%% Xor Mask functions

% this function creates a xor mask 
function seq = GenGoldCode(Nc, L, seed)
arguments
    Nc (1,1) = 1600
    L (1,1) = 7200
    seed (1,1) = 0x12345678
end

x1 = zeros(Nc + L + 31, 1);
x2 = zeros(Nc + L + 31, 1);
x2(1:32) = flip(fdec2bin(seed, 32));
x1(1) = 1;
for n = 1:(Nc + L)
    x1(n + 31) = xor(x1(n + 3), x1(n));
    x2(n + 31) = xor(xor(x2(n + 3), x2(n + 2)), xor(x2(n+1), x2(n)));
end
seq = xor(x1(Nc + 1:Nc+L), x2(Nc + 1:Nc+L));

end

% this function converts decimal to bits
function decimal=fdec2bin(dec_num, min_dig)
    decimal=dec2bin(dec_num,min_dig)-'0';
end

%%  Cyclic Buffer functions

% this function removes a cyclic buffer
function BitsAfterCycleBufferRmoval = CyclicBufferRemover(bits,start,buffer_size)
% bits - bits of info
% start - start of info
% buffer_size - size of info

BitsAfterCycleBufferRmoval=[bits(start:buffer_size),bits(1:start-1)];

end

%% Turbo Decoder functions
function [DeInterlacedP1, DeInterlacedP2] = Deinterlacer(InterleavedInerlacedP)
DeInterlacedP1= InterleavedInerlacedP(1:2:end);
DeInterlacedP2= InterleavedInerlacedP(2:2:end);
end

function DummySpots=Dummy(SubBlockLength, DummySize, DummyMatSize, DummyFormat)

RandomMat=(1:SubBlockLength);

%insert Dummies
DummyBlock=[NaN(1,DummySize),RandomMat];

%reshape to mat by rows
DummyBlock=reshape(DummyBlock,flip(DummyMatSize)).';

%mix the columns
DummyBlock=DummyBlock(:,DummyFormat);

%reshape to row by columns
DummyBlock=reshape(DummyBlock,1,[]);

%remove Dummies
DummySpots=find(isnan(DummyBlock));

end

function DeDummyedBlock=DeDummy(InterleavedSubBlock, DummySize, DummyMatSize, DummyFormat, DummySpots)

%insert Dummies
InterleavedSubBlockLen = length(InterleavedSubBlock) + DummySize;
DeDummyedBlock = nan(1, InterleavedSubBlockLen); 
data_idx = setdiff(1:InterleavedSubBlockLen, DummySpots); %find the indexes that are not nan in the dummyd data
DeDummyedBlock(data_idx) = InterleavedSubBlock; 

%reshape to mat by columns 
DeDummyedBlock=reshape(DeDummyedBlock,DummyMatSize);

%mix back the columns
DeDummyedBlock=DeDummyedBlock(:,DummyFormat);

%reshape to row by rows
DeDummyedBlock=reshape(DeDummyedBlock.',1,[]);

%remove Dummies and cut padding
DeDummyedBlock=DeDummyedBlock(DummySize+1:end-4);

end

function DecodedData=TurboDecoder(S, P1, P2 ,interliver, ConvCodeRate)

ConstraintLen= 4; % because the polynomial degree is 3 
FeedBackConn =13; % because 13 in octal base is 11 which is 1011 in binary base which is equal to 1+D^2+D^3 which is the polynom of the encryption
CodeGen= [13 15]; % 13[dec]= 1101[bin]=1+D+D^3, 15[dec]=1111[bin]=1+D+D^2+D^3 (G)
trellis = poly2trellis(ConstraintLen,CodeGen,FeedBackConn); 

tbdepthPres=[1/2 5;2/3 7.5;3/4 10;5/6 15]; % common tbdepth by CodeRate according to matlab documentation
tbdepth=tbdepthPres(tbdepthPres(:,1)==ConvCodeRate,2)*(ConstraintLen-1); %because rate 1/2 code has a traceback depth of 5 × (ConstraintLength – 1)
opmode='trunc'; % because the encoder is assumed to have started at the all-zeros state.
dectype='hard'; % because the input values are 0 or 1.

S= vitdec(Interlacer(S, P1) ,trellis,tbdepth,opmode,dectype); % interlace S&P1 and operate viterbi on the combination
S=Interliver(S, interliver);
DecodedData= vitdec(Interlacer(S, P2),trellis,tbdepth,opmode,dectype); % interlace S after viterbi with P1 and interliving & P2 and operate viterbi on the combination
DecodedData=DeInterliver(DecodedData, interliver);

end

function InterlacedBlocks = Interlacer(B1, B2)
%sort of an interliver ,used to interlace P and S 
B1 = B1(:).';
B2 = B2(:).';
BlockMat = [B1 ; B2];
InterlacedBlocks = BlockMat(:).';
end

function InterlivedBlock = Interliver(UnInterlivedBlock, interleaver)
InterlivedBlock = UnInterlivedBlock(interleaver + 1);
end

function DeInterlivedBlock = DeInterliver(InterlivedBlock, interliver)
interliver_mapping = Interliver(1:length(interliver), interliver);
[~, reverse_mapping] = sort(interliver_mapping);
DeInterlivedBlock = InterlivedBlock(reverse_mapping);
end


%% CRC functions

function IsValidCRC = CRC_Vailidation(deECCbits, deECCbitsLen, polynom,init)

% validates the check value of the Cyclic redundancy check 
% parameters:
% deECCbits - the bits after decoding the xor mask, cyclic buffer and error correction code (interliver)
% deECCbitsLen - the length of the bit sequence (data+CRC_CheckValue)
% polynom - the characteristic polynomial
% init - initialize the crc value before the polynom division

residue = zeros(1,length(init));

for i=0:deECCbitsLen - 1
    MaxIdx = residue(end);
    residue = [deECCbits(end-i),residue(1:end-1)];
    if MaxIdx == 1
        residue =  xor(residue,polynom(1:end-1)); % modulate with polynom
    end
end

IsValidCRC = (residue== init);

end

%% parsing function
function parsing(deECCbits)