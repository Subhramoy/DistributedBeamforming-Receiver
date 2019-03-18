clear all; close all; clc;  %#ok

%% Basic Parameters
SNR_desired = 30;
Atx = 3;
Arx = 3;
GainTx = 30;
GainRx = 30;
Ampl = 0;
Noise = 95;

%% Compute Power limits in Base Band
load('centralized/data/information4.mat','payload1');
w_max_real = 1/max(real(payload1));
w_max_imag = 1/max(real(payload1));
w_max = w_max_real + 1i*w_max_imag;
BBPowMax = w_max*w_max';  % I and Q values cannot exceed 1 in Base Band (DAC limitation)
BBPowPayload = (payload1'*payload1)/length(payload1);

d_max = f_BFLimits1(SNR_desired,Atx,Arx,GainTx,GainRx,Ampl,Noise,BBPowMax,BBPowPayload)

function d_max = f_BFLimits1(SNR_desired,Atx,Arx,GainTx,GainRx,Ampl,Noise,BBPowMax,BBPowPayload)

    N = 4;  % Number of antennas used in Beamforming
    f = 900e6;  % in Hertz
    K = 20*log10(f) - 147.55;

    P_tx_max = pow2db(BBPowPayload * BBPowMax) + 30;  % in dBm
    FSPL_max = -SNR_desired + (P_tx_max + 20*log10(N) + Atx + Arx + GainTx + GainRx + Ampl - Noise);
%     FSPL_max = -FSPL_max;

%     FSPL_max = 20log10(d_max) + K
    d_max = 10^((FSPL_max + K)/20);  % in meters

end