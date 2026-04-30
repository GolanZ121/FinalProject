%% SDR side
clc; clear; close all;

%% Params
NFFT = 1024;
scs = 15e3;
CorrectFs = scs * NFFT;
Fs = 15.36e6;
[p,q] = rat(CorrectFs / Fs);
cps_lens = [80 72 72 72 72 72 72 72 80];
ZC1_sym = 4;
ZC2_sym = 6;
SPP = sum(cps_lens) + length(cps_lens) * NFFT;  % Samples per packet
Fcs = [2.399 2.4145 2.4259 2.4445 2.4595] * 1e9;
FcIdx =1;

start_sc = 213;
stop_sc = 813;
%% create ZadoffChu sequences
zc1 = create_zc_OFDM(600);
zc2 = create_zc_OFDM(147);

%%
rx = comm.SDRuReceiver( ...
    'Platform', 'B210', ...
    'SerialNum', '34D62B0', ...
    'CenterFrequency',    2.4295e+09, ...
    'Gain', 30, ...
    'MasterClockRate', Fs, ...
    'SamplesPerFrame', 64000, ...
    'DecimationFactor', 1);

% subplot(2,1,1);
%%

num_frames = 50;

while true
    data_vec = [];
    for i = 1:num_frames
        data = double(rx());
        data_vec = [data_vec; data(:)];
    end
    data_res = resample(data_vec, p, q);
    
    % subplot(2,1,1);
    [time_c_rough, time_idx_rough] = cross_corr(zc1, data_res);
    threshold = 40* mean(abs(time_c_rough));
    threshold = 50e3;
    samples_before_zc1 = sum(cps_lens(1:ZC1_sym)) + (ZC1_sym-1)*NFFT;
    start_index = time_idx_rough - samples_before_zc1;
    plot(abs(time_c_rough));
    yline(threshold, "--r");
    if max(abs(time_c_rough)) > threshold
        break
    else
        FcIdx = FcIdx + 1;
        FcIdx = mod(FcIdx, length(Fcs)) + 1
        rx.CenterFrequency = Fcs(FcIdx);
    end
    % subplot(2,1,2);
    % spectrogram(data_res,100,80,100,Fs, 'centered', 'yaxis' );
    drawnow;
end

%%
figure;
spectrogram(data_res,100,80,100,Fs, 'centered', 'yaxis' );
%% Time Sync (Rough because we still have the freq shift)
[time_c_rough, time_idx_rough] = cross_corr(zc1, data_res);
samples_before_zc1 = sum(cps_lens(1:ZC1_sym)) + (ZC1_sym-1)*NFFT;
start_index = time_idx_rough - samples_before_zc1;
figure;
plot(abs(time_c_rough));

%% Cut packet
figure;
time_aligned_packet = data_res(start_index:start_index + SPP - 1);
spectrogram(time_aligned_packet,100,80,100,Fs,'centered', 'yaxis')

%% Coarse frequency correction (using CPs)
coarse_freq_offset = find_freq_offset(time_aligned_packet, cps_lens, NFFT, Fs);

% === > Fix the freq shift < === %
data_vec = fix_freq(data_vec, coarse_freq_offset, Fs);

%% Time Sync (Rough because we still have the freq shift)
[time_c_rough, time_idx_rough] = cross_corr(zc2, data_vec);
samples_before_zc1 = sum(cps_lens(1:ZC2_sym)) + (ZC2_sym-1)*NFFT;
start_index = time_idx_rough - samples_before_zc1;
figure;
plot(abs(time_c_rough));

%% Cut packet
figure;
time_aligned_packet = data_vec(start_index:start_index + SPP - 1);
spectrogram(time_aligned_packet,100,80,100,Fs,'centered', 'yaxis')

%% Estimate channel
recievedZC = time_aligned_packet(NFFT*(ZC1_sym-1) + sum(cps_lens(1:4)) + 1: NFFT*ZC1_sym + sum(cps_lens(1:4)));
H = fft(recievedZC) ./ fft(zc1);


demodulated_data = ofdm_demod(time_aligned_packet,cps_lens(2),cps_lens(1),NFFT,start_sc,stop_sc, length(cps_lens), H);



%% Functions

function fixed_data = fix_freq(data, f_offset, Fs)
CORRECT_time_axis = (0:length(data)-1).' / Fs;
fixed_data = data .* exp(-1j * 2 * pi * f_offset * CORRECT_time_axis);
end

function freq_offset = find_freq_offset(packet, cps_lens, NFFT, Fs)
cp_start_index = 1;
num_of_syms = length(cps_lens);
vec = zeros(num_of_syms, 1);
for i = 1: num_of_syms
    cp_i1 = double(packet(cp_start_index : cp_start_index + cps_lens(i) - 1));
    cp_i2 = double(packet(cp_start_index + NFFT : cp_start_index + NFFT + cps_lens(i) -1));
    vec(i) = sum(cp_i1 .* conj(cp_i2));
    cp_start_index = cp_start_index + NFFT + cps_lens(i);
end

freq_offset = -angle(mean(vec)) * Fs / (2 * pi * NFFT);
end

function [c, idx] = cross_corr(sig1, sig2)
c = conv(sig2, flip(conj(sig1)), "valid");
[~, idx] = max(abs(c));
end

function ZCOFDM = create_zc_OFDM(root)
ZC = zadoffChuSeq(root, 601);
ZCPAD = [zeros(212, 1); ZC; zeros(211,1)];
ZCOFDM = ifft(ifftshift(ZCPAD));
end

function data = load_samples(filename, file_Fs, Fs)
fid = fopen(filename , "r");
samples = fread(fid, [2 inf], "float32=>double");
fclose(fid);
% data = complex(samples(1,:), samples(2, :));
data = resample(complex(samples(1,:), samples(2, :)).', Fs, file_Fs);
end
function [X, freqs] = plot_fft(signal, Fs)
arguments
    signal
    Fs
end
X = fftshift(fft(signal));
freqs = linspace(-Fs/2, Fs/2, length(X));
plot(freqs, (abs(X)));
title("FFT");
xlabel("Frequency [Hz]");
ylabel("Magnitude [-]");
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