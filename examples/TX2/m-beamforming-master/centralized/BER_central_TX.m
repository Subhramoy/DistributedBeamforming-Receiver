%% Clean environment
if exist('transmitter','var')
    release(transmitter);
end
clear all; close all; clc;  %#ok
addpath('..','-end');  % add shared code between centralized and decentralized

%% Configure experiment
numTxAntennas = 2;  % Select between 1, 2 and 4
maxIter       = 30000;  % Maximum transmissions over the air
gain          = 8; % Default is 8 (short distance - indoor)
fileName      = 'weights.bin';  % File location for channel estimation (decentralized)
verbose       = true;  % debugging
bfAlg         = 'RWA';  % options are: basic, RWA
%% Comms parameters
Atx = 0;  % Antenna Gain at the transmitter (in dB)
Noise = 45;  % Noise floor (in dBm)
SINRdb = 20;  % SNR target at the receiver (in dB)
%% Water-filling configuration
DeltaP = 0;  % robustness factor (power in linear)
inc = 0.001;  % water filling increment

%% Configure radios
fprintf('LOG - Configuring experiments...\n');
[radioConfig,bits,~,payload1,trainingSig,fid,channelEst_old] = ...
                                   f_loadConfig_TX(numTxAntennas,fileName);


%% Create Transmitter
fprintf('LOG - Creating TX System Object. Please wait...\n');
transmitter = comm.SDRuTransmitter('Platform',radioConfig.txPlatform,...
                                   'IPAddress',radioConfig.txIPAddrs,...
                                   'MasterClockRate',radioConfig.txMasterClockRate,...
                                   'InterpolationFactor',radioConfig.txInterpolationfactor,...
                                   'ChannelMapping',radioConfig.ChannelMapping,...
                                   'CenterFrequency',900e6,...
                                   'Gain',gain,...
                                   'ClockSource','External',...
                                   'PPSSource','External');

%% Compute power limitations for Waterfilling
w_max_real = 1/max(real(payload1));
w_max_imag = 1/max(real(payload1));
w_max = w_max_real + 1i*w_max_imag;
BBPowMax = w_max*w_max';  % I and Q values cannot exceed 1 in Base Band (DAC limitation)
BBPowPayload = (payload1'*payload1)/length(payload1);
                               
%% Main loop
fprintf('LOG - Begin transmission over main loop...\n');
for i = 1:maxIter
    fseek(fid,0,'bof'); % Read from the beginning of the file
    channelEst = fread(fid,numTxAntennas*2,'double');  % Complex = 2 double

    if length(channelEst) == numTxAntennas*2
        % File content has expected length
        channelEst = (   channelEst(1:numTxAntennas) + ...
                       1j*channelEst(numTxAntennas+1:numTxAntennas*2)).';
        channelEst_old = channelEst;
    else
        % Use last feedback
        channelEst = channelEst_old;
    end

    % Compute beamforming weights
    if strcmp(bfAlg,'basic') % Channel inversion
        beamWeight = f_BFPayload(channelEst);
    elseif strcmp(bfAlg,'RWA')  % Channel inversion + Waterfilling
        [beamWeight,~] = f_BFPayload_WF(channelEst,DeltaP,inc,SINRdb,0,0,Noise,BBPowMax,BBPowPayload);
    end

    % Apply beamforming to payload
    payload = payload1*beamWeight;
    
    if verbose
        fprintf('Iter %d - ',i);
        for id = 1:numTxAntennas
            fprintf('H(%d) = %.7f + %.7fj | ',id,real(channelEst(id)),imag(channelEst(id)));
        end
        fprintf('\n');
        for id = 1:numTxAntennas
            fprintf('W(%d) = %.7f + %.7fj | ',id,real(beamWeight(id)),imag(beamWeight(id)));
        end
        fprintf('\n');
    end

    % Send signals to the radios
    txSig = [trainingSig; zeros(400,numTxAntennas); payload; zeros(500,numTxAntennas)] * 0.2;
    transmitter(txSig);
end

release(transmitter);  % Release the System Object for future use
fclose(fid);  % close file descriptor for channel feedback


% EOF