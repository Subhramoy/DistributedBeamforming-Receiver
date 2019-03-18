  %% Clean environment
if exist('transmitter','var')
    release(receiver);
end
clear all; close all; clc;  %#ok

%% Configuration
maxIter      = 300;  % Number of iterations

%% Configure radios
fprintf('LOG - Configuring experiments...\n');
[radioConfig,bits,~,~,trainingSig,fid,channelEst_old] = ...
                                   f_loadConfig_RX(1,'rubish');

% Connect to radio
fprintf('LOG - Creating RX System Object. Please wait...\n');
receiver = comm.SDRuReceiver('Platform',radioConfig.rxPlatform1,...
                             radioConfig.rxIDProp1,radioConfig.rxID1,...
                             'MasterClockRate',radioConfig.rxMasterClockRate1,...
                             'DecimationFactor',radioConfig.rxDecimationfactor1,...
                             'ClockSource','External',...
                             'Gain',0,...
                             'CenterFrequency',900e6,...
                             'EnableBurstMode',true,...
                             'SamplesPerFrame',50e3,...
                             'OutputDataType','double',...
                             'NumFramesInBurst',0);
                         
%% Main loop
fprintf('LOG - Begin reception over main loop...\n');
noise_power_lin = zeros(1,maxIter);
for i = 1:maxIter
    [rxSig, len] = receiver();
    noise_power_lin(i) = std(rxSig)^2;
    fprintf('Noise Power = %.5f (dBm)\n',pow2db(noise_power_lin(i)));
end
noise_power_av_lin = mean(noise_power_lin);
noise_power = pow2db(noise_power_av_lin) + 30;
fprintf('   Noise Power = %.5f (dB)\n',pow2db(noise_power_av_lin));