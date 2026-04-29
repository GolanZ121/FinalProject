%% sync_samples.m
% Returning the fully time and freq synced packet

clc; clear; close all;

%% Params
filename = 'rec1.mat';
file_Fs = 15.36e6;
NFFT = 1024;
scs = 15e3;
Fs = scs * NFFT;
cps_lens = [80 72 72 72 72 72 72 72 80];
ZC1_sym = 4;
ZC2_sym = 6;
SPP = sum(cps_lens) + length(cps_lens) * NFFT;  % Samples per packet

%% Load samples 
data = load_samples(filename, file_Fs, Fs);
data(isnan(data)) = 0;

load(filename);
data = data2;

%% create ZadoffChu sequences
zc1 = create_zc_OFDM(600);
zc2 = create_zc_OFDM(147);

%% Time Sync (Rough because we still have the freq shift)
[time_c_rough, time_idx_rough] = cross_corr(zc1, data);
samples_before_zc1 = sum(cps_lens(1:ZC1_sym)) + (ZC1_sym-1)*NFFT;
start_index = time_idx_rough - samples_before_zc1;
plot(abs(time_c_rough));
time_aligned_packet = data(start_index:start_index + SPP - 1);

%% Coarse frequency correction (using CPs)
coarse_freq_offset = find_freq_offset(time_aligned_packet, cps_lens, NFFT, Fs);

% === > Fix the freq shift < === %
data = fix_freq(data, coarse_freq_offset, Fs);


%% Fine time sync 
[time_c_fine, time_idx_fine] = cross_corr(zc1, data);
start_index = time_idx_fine - samples_before_zc1;

time_aligned_packet = data(start_index:start_index + SPP - 1);

%% Estimate channel
recievedZC = time_aligned_packet(NFFT*(ZC1_sym-1) + sum(cps_lens(1:4)) + 1: NFFT*ZC1_sym + sum(cps_lens(1:4)));
H = fft(recievedZC) ./ fft(zc1);

%% Save Synced samples and channel as files
fid = fopen('output.32fc','wb');       
x = single(time_aligned_packet); 
interleaved = zeros(2*numel(x),1,'single');
interleaved(1:2:end) = real(x); 
interleaved(2:2:end) = imag(x); 
fwrite(fid, interleaved, 'float32');
fclose(fid);

save("channel.mat", "H");

%% 

function fixed_data = fix_freq(data, f_offset, Fs)
    CORRECT_time_axis = (0:length(data)-1).' / Fs;
    fixed_data = data .* exp(-1j * 2 * pi * f_offset * CORRECT_time_axis);
end

function freq_offset = find_freq_offset(packet, cps_lens, NFFT, Fs)
    cp_start_index = 1;
    num_of_syms = length(cps_lens);
    vec = zeros(num_of_syms, 1);
    for i = 1: num_of_syms
        cp_i1 = packet(cp_start_index : cp_start_index + cps_lens(i) - 1);
        cp_i2 = packet(cp_start_index + NFFT : cp_start_index + NFFT + cps_lens(i) -1);
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