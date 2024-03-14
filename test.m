clearvars;
close all;

seed = rng(51342);
N=32;
M=128;
lmo = 3;
countarr = zeros(3, 9);
Pav = 1;
count = 0;
errorcal = comm.ErrorRate;
%% channel options
EVA = [[0 30 150 310 370 710 1090 1730 2510]*1e-9;
    [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9]];
ETU = [[0 50 120 200 230 500 1600 2300 5000]*1e-9;
    [-1.0 -1.0 -1.0 0.0 0.0 0.0 -3.0 -5.0 -7.0]];
% [~, ~, ~, channelLen] = OTFS.gen_chan_param(M, N, 15e3, 1/15e3, 5.9e9, 500, EVA);
[~, ~, ~, channelLen] = OTFS.gen_chan_param(M, N, 15e3, 1/15e3, 5.9e9, 500, ETU);
L = channelLen + 1 * lmo;
if ~mod(L, 2)
    L = L + 1;
end
%% equal seq length
L = 15;

%% generate 3 frames adn channel
R = zeros(M, 3*N);
noisePower = Pav / 10^(10/10);
CFO = rand - 0.5; % (-0.5~0.5)
[delay, doppler, ch] = OTFS.gen_chan_param(M, N, 15e3, 1/15e3, 5.9e9, 500, EVA);
% [delay, doppler, ch] = OTFS.gen_chan_param(M, N, 15e3, 1/15e3, 5.9e9, 500, ETU);
intDoppler = round(doppler);
H = OTFS.gen_channel(M, N, delay, intDoppler, ch, 0, false);
noise = sqrt(noisePower/2) * (randn(size(1, M*N)) + 1i*randn(size(1, M*N)));

for three = 0:2
    [~, moddata] = OTFS.gen_bitstream(M, N, L, 4);
    [xdd, pilotseq, mp] = OTFS.gen_data_grid(N, 'pcp', noisePower, L, moddata, 10);    
    s = ifft(xdd.')*sqrt(N);
    s = s(:);
    r = H*s + noise;
    r_prime = r .* exp(1i*2*pi*CFO*(1:M*N)/M/N);
    R(:, 1+three*N:(three+1)*N) = reshape(r_prime, M, N);
end

%% p_d[m]
