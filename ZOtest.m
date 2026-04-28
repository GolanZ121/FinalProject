clc;
clear; 
close all;
NFFT = 1024;
og_Fs = 15.36e6;



%% Generate ZadoffChu
ZC600 = zadoffChuSeq(600, 601);
ZC600PAD = [zeros(212, 1); ZC600; zeros(211, 1)];
ZC600OFDM = ifft(ifftshift(ZC600PAD));

f_shift = 5e3;

Ts = (NFFT + 80)/ og_Fs;
time = 0: 1/og_Fs: Ts - 1/og_Fs;
KNOWNCP = ZC600OFDM(NFFT-80 + 1: NFFT);

ZC_with_cp = [ZC600OFDM(NFFT-80 + 1: NFFT); ZC600OFDM];
shiftedZO = ZC_with_cp .* exp(1j * 2 * pi * f_shift * time).';
tZO = [zeros(10000, 1); shiftedZO; zeros(234235, 1)];
Ttotal = length(tZO) / og_Fs;
%%
figure;
spectrogram(tZO,64,32,1024,og_Fs,'yaxis', 'centered')

[c1, idx1] = cross_corr(ZC600OFDM, tZO);
figure;
plot(1:length(c1), [abs(c1)]);

% Get the ZO base on the time sync
ZOcut = tZO(idx1-80:idx1 + NFFT-1);
% Fix the Coarse frequency shift
time = 0: 1/og_Fs: Ttotal - 1/og_Fs;
%% Find coarse frequency offset

cp11 = ZOcut(1:80);
cp12 = ZOcut(NFFT+1: NFFT + 80);
freq_offset = mean(angle(cp11 .* conj(cp12)))  * og_Fs / (2*pi*(NFFT + 80));

tZO = tZO .* exp(1j * 2 * pi * freq_offset * time).';

%%
[c1, idx1] = cross_corr(ZC600OFDM, tZO);
figure;
plot(1:length(c1), [abs(c1)]);

% Get the ZO base on the time sync
ZOcut = tZO(idx1-80:idx1 + NFFT-1);

%% functions

function [conv_res, max_idx] = cross_corr(sig2, sig1)
    conv_res = conv(sig1, flip(conj(sig2)), "valid");
    [~, max_idx] = max(conv_res);
end

function [X, freqs] = plot_fft(signal, Fs, nfft)
arguments
    signal 
    Fs 
    nfft = 1024
end
    X = fftshift(fft(signal, nfft));
    freqs = linspace(-Fs/2, Fs/2, length(X));
    figure;
    plot(freqs, (abs(X)));
end


function [X, freqs] = zoomed_ft(x, f_start, f_end, Fs, res)
    NFFT = Fs / res;

    %% =========== Z TRANSFTORM =========== %%
        m = NFFT;
        w = exp(-1j*2*pi*(f_end-f_start)/(m*Fs)); % Spiral ratio
        a = exp(1j*2*pi*f_start/Fs); % Starting point
        
        % Compute the zoom spectrum
        X = czt(x, m, w, a);
        freqs = linspace(f_start, f_end, m);

end