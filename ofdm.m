%% 
%% clean
clc; clear; close all;

%%
load("output.mat");
load("channel.mat");
%% Scripts to fix the phase 

% Params
Fs = 15.36e6;
Ncp = 72;
ex_Ncp = 80;
Nsc = 1024;
start_sc = 213;
stop_sc = 813;
num_of_symbols = 9;
NFFT = 1024;


%% Script to examine ofdm function

demodulated_data = ofdm_demod(data,Ncp,ex_Ncp,Nsc,start_sc,stop_sc, num_of_symbols, H);


%% functions

% this function searches the phase diff
function phase = phase_sync(block4,block6)
    %block4 - first Zc symbol
    %block6 - second Zc symbol

    % create the 2 ZadoffChu symbols
    ZC600 = zadoffChuSeq(600, 601);
    ZC600PAD = [zeros(212, 1); ZC600; zeros(211, 1)];
    ZC600OFDM = ifft(ifftshift(ZC600PAD));
    
    ZC147 = zadoffChuSeq(147, 601);
    ZC147PAD = [zeros(212, 1); ZC147; zeros(211, 1)];
    ZC147OFDM = ifft(ifftshift(ZC147PAD));

    % calculate the phase diff
    phase1 = angle(conj(block4)*ZC600OFDM);
    phase2 = angle(conj(block6)*ZC147OFDM);
    phase = (phase1 + phase2)/2;
end

% this function demodulates ofdm
function data = ofdm_demod(packet,Ncp,ex_Ncp,Nsc,start_sc,stop_sc, num_of_symbols, H)
    % packet - the packet of information
    % Ncp - cyclic prefix, a partial copy of the end of the symbol to create  a
    % correlation peak at each start of symbol
    % ex_Ncp - extended cp
    % start_sc - is the start of the active subcarrier range
    % stop_sc - is the start of the active subcarrier range
    % num_of_symbols - a given information about num of symbols in a packet
    
    clean_packet = remove_cp(packet,Ncp,ex_Ncp,Nsc,num_of_symbols); % packet without cp
    
    slot = reshape(clean_packet,Nsc,[]); % create the slot matrix
    slot_fft = fftshift(fft(slot) ./ repmat(H, 1, 9) ,1)./(sqrt(Nsc)); % according to the formula (not sure if needed in our case)
    slot_fft = slot_fft(start_sc:stop_sc,:); % cut the irrelevant frequencies (frequency domain)
    slot_fft = [slot_fft(1:300,:);slot_fft(302:end,:)]; % the mid subcarrier is null
    slot_fft = [slot_fft(:,2:3),slot_fft(:,5),slot_fft(:,7:end)]; % The data is stored in symbols [2,3,5,7,8,9]
    
    figure;
    title("after Channel Est")
    plot(slot_fft(:, :), "o");
    
    % transform the slot to bits of data
    data = slot_fft(:);
    data = pskdemod(data,4,0.25*pi);
    data(data == 2) = 5;
    data(data == 1) = 2;
    data(data == 5) = 1;
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
    clean_packet = zeros(1,Nsc*num_of_symbols);
    symbol_index = Ncp + 1;
    % at each iteration of the loop I cut the cp and paste the data
    for i = 1:num_of_symbols
        %symbols 1 and 9  have extended cp
        if i == 1 || i == num_of_symbols
            symbol_index = symbol_index + cp_diffs ;
        end
        end_symbol_index = symbol_index + Nsc  - 1;
       clean_packet((i-1)*Nsc + 1:i*Nsc) = packet(symbol_index :end_symbol_index ) ;
       symbol_index = symbol_index + Nsc + Ncp; % calculate the next symbol index
    end
end