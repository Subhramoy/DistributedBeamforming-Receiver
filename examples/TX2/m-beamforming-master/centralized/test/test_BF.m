clear all; close all; clc;
addpath('../data_experiments','-end');  % to Load real channels
addpath('../data','-end');  % to Load Payload
addpath('../utils','-end');  % load utility functions
addpath('../../','-end');  % load main functions (Beamforming)

%% Comms parameters
Atx = 3;  % in dB
Gtx = 20;  % in dB
Noise = 105;  % in dBm
SNR_demand = 25;  % in dB
P_ch_desired = 120;  % emulate channel losses
%% Water-filling configuration
DeltaP = 0;  % robustness factor (power in linear)
inc = 0.001;  % water filling increment
update = 20;

%% Load channel
load('expOutdoor_4Tx1Rx_TxHov_RxStat_-20dBm_1-5m.mat','chTot');
N_ant = size(chTot,1);
nanLocations = find(isnan(chTot(1,:)));
chTot(:,nanLocations) = [];
Niters = size(chTot,2);

%% Load Payload
load('information4.mat','payload1');
BBPowPayload = (payload1'*payload1)/length(payload1);

%% Compute demanded Tx power in Base Band (BB)
P_tx_demand = SNR_demand + P_ch_desired - Atx - Gtx - Noise;  % in dBm
fprintf('LOG - Requested Ptx: %.3f (dBm)\n',P_tx_demand);
fprintf('LOG - Requested SNR: %.3f (dB)\n',SNR_demand);

%% Compute maximum Tx power offered in Base Band (BB)
w_max = 2 + 2j;  % given current txSig
BBPowMax = BBPowPayload*(w_max*w_max');
BBPowMax_rep = repmat(BBPowMax,1,N_ant);
P_tx_offered = pow2db(sum(BBPowMax_rep)) + 30;  % in dBm
fprintf('LOG - Offered MAX Ptx: %.3f (dBm)\n',P_tx_offered);
SNR_offered = P_tx_offered + Atx + Gtx - P_ch_desired + Noise;  % in dB
fprintf('LOG - Offered MAX SNR: %.3f (dB)\n',SNR_offered);

%% Generate output variables
P_tx_basic = zeros(Niters,1);
SNR_basic = zeros(Niters,1);
P_tx_RWA = zeros(Niters,1);
SNR_RWA = zeros(Niters,1);
P_tx_RWA_R = zeros(Niters,1);
SNR_RWA_R = zeros(Niters,1);
for iter = 1:Niters
    %% Retrieve channel information
    ch = chTot(:,iter).';  % pre-store channel in [1 x M] dimensions (compatibility)
    P_ch = -pow2db(ch*ch');  % compute channel power
    ch = sqrt(db2pow(-P_ch_desired)/db2pow(-P_ch)).*ch;
    P_ch = -pow2db(ch*ch');  % compute NEW channel power

    %% Basic beamforming
    w_basic = f_BFPayload(ch);
    P_tx_basic(iter) = pow2db(BBPowPayload*(w_basic*w_basic')) + 30;  % in dBm
    SNR_basic(iter) = P_tx_basic(iter) + Atx + Gtx - P_ch + Noise;  % in dB

    %% Waterfilling beamforming
    [w_RWA,~] = f_BFPayload_WF(ch,DeltaP,inc,SNR_demand,Atx,Gtx,Noise,BBPowMax,BBPowPayload);
    P_tx_RWA(iter) = pow2db(BBPowPayload*(w_RWA*w_RWA')) + 30;  % in dBm
    SNR_RWA(iter) = P_tx_RWA(iter) + Atx + Gtx - P_ch + Noise;  % in dB
    
    %% Waterfilling beamforming - Robust
    if mod(iter,20)==1
        [w_RWA_R,~] = f_BFPayload_WF(ch,DeltaP,inc,SNR_demand,Atx,Gtx,Noise,BBPowMax,BBPowPayload);
        w_RWA_abs_old = abs(w_RWA_R);
    else
        w_temp = f_BFPayload(ch);
        angles = angle(w_temp);
        w_RWA_R = w_RWA_abs_old.*(cos(angles)+...  % Complex beamWeights
                                  1i.*sin(angles));
    end
    P_tx_RWA_R(iter) = pow2db(BBPowPayload*(w_RWA_R*w_RWA_R')) + 30;  % in dBm
    SNR_RWA_R(iter) = P_tx_RWA(iter) + Atx + Gtx - P_ch + Noise;  % in dB
    
    %% Plotting - sample for weight distribution
    if iter==1
        NumStackElements = 2;  % basic and RWA
        NumGroupsPerAxis = 4;  % number of channels
        NumStacksPerGroup = 2;  % channel and weights
        groupLabels = {'h_1';'h_2';'h_3';'h_4'};
        A(:,1,1) = 1./(ch/max(abs(ch))).';
        A(:,1,2) = w_basic;
        A(:,2,1) = A(:,1,1);
        A(:,2,2) = w_RWA;
        plotBarStackGroups(abs(A), groupLabels, 1);
    end
end

%% Logging
fprintf('LOG - Offered Basic Ptx: %.3f (dBm)\n',mean(P_tx_basic));
fprintf('LOG - Offered Basic SNR: %.3f (dB)\n',mean(SNR_basic));
fprintf('LOG - Offered Waterfilling Ptx: %.3f (dBm)\n',mean(P_tx_RWA));
fprintf('LOG - Offered Waterfilling SNR: %.3f (dB)\n',mean(SNR_RWA));

%% Plotting - average
figure; hold on;
plot(1:Niters,P_tx_basic,'color','b');
plot(1:Niters,P_tx_RWA,'color','r');
xlabel('Iterations');
ylabel('Power (dBm)');
title('Transmitted power');
legend('Basic BF','Waterfilling BF');
pos = get(gcf, 'Position');
set(gcf,'Position',[pos(1) pos(2) 347 159]);
figure; hold on;
plot(1:Niters,SNR_basic,'color','b');
plot(1:Niters,SNR_RWA,'color','r');
xlabel('Iterations');
ylabel('SNR (dB)');
title('Achieved SNR');
legend('Basic BF','Waterfilling BF');
pos = get(gcf, 'Position');
set(gcf,'Position',[pos(1) pos(2) 347 159]);