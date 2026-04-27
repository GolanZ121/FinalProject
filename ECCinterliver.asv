clc; clear; close all;

%% definitions
Fs = 10e6;

% load the raw_samples file 
filename = 'C:\Users\ZBOOK\FinalProject\processed_bits\processed_bits12_Feb_2026_09_30_39_442.txt';
fid = fopen(filename, 'r');

% Read all data
ProccesedBits = dec2bin(fread(fid)); %len=4236
fclose(fid);

% Turbo Decoder parameters
interliver= mod((43 * (1/1408) + 88 * (1/1408).^2), 1408);
[q,r] = polydiv([1 1 1 1],[1 0 1 1]);
G=[eye, q];
ConvCodeRate=1/2;
FinalCodeRate=1/3;
SubBlockLength=1412;
DummyFormat=[1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31,2,18,10,26,6,22,14,30,4,20,12,28,8,24,16,32];
DummySize=28;
DummyMatSize=[45, 32];

% test
data=TurboDecoder(ProccesedBits, SubBlockLength, interliver, DummySize, DummyMatSize, DummyFormat);

%% Turbo Decoder

function data=TurboDecoder(ProccesedBits, SubBlockLength, interliver, DummySize, DummyMatSize, DummyFormat)
    %input: parameters and processed (after xor and cyclic shift removal) bits
    %output: data bits for CRC 
    [S, P1, P2]= PreDecoderDecode(ProccesedBits, SubBlockLength, interliver, DummySize, DummyMatSize, DummyFormat);
    
    %viterbi parameters:
    num_states=2^(1408-1);
    states=(1:num_states);
    % InitialProbabilities=ones(sqrt(num_states))/2;
    % for s = 1:num_states
    %     ObservationsProbabilitiesMat(s, 1) = InitialProbabilities(s) * ObservationsProbabilitiesMat(s, ObservationsSequnce(1));
    % end
    % TransitionProbabilitiesMat = 
    % ObservationsSequnce =

    init=;
    trans=;
    emit=;
    obs=;

    data=Viterbi(states,init,trans,emit,obs);

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

function path=Viterbi(States,InitialStateProbabilities,TransitionStatesMatrix,EmissionMatrix,ObservationsSequnce)
    
    NumOfStates=length(States);
    ObservationSequnceSize=length(ObservationsSequnce);

    probabilities=zeros(ObservationSequnceSize,NumOfStates);
    previousStateMat=NaN(ObservationSequnceSize,NumOfStates);

    for s=1:NumOfStates
        probabilities(0,s)=InitialStateProbabilities(s)*EmissionMatrix(s,ObservationsSequnce(0));
    end

    for t=1:ObservationSequnceSize-1
        for s=1:NumOfStates
            for r=NumOfStates
                new_prob=probabilities(t-1,r)*TransitionStatesMatrix(r,s)*EmissionMatrix(s,obt(t));
                if new_prob>probabilities(t,s)
                    probabilities(t,s)=new_prob;
                    previousStateMat(t,s)=r;
                end
            end
        end
    end

    path=NaN(ObservationSequnceSize,1);
    m=max(probabilities(ObservationSequnceSize,States));
    path(end)=States(probabilities(ObservationSequnceSize,States)==m);

    for t=ObservationSequnceSize-2:0
        path(t)=previousStateMat(t+1,path(t+1));
    end
end

