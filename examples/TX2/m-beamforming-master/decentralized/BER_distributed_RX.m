  %% Clean environment
if exist('transmitter','var')
    release(receiver);
end
clear all; close all; clc;  %#ok
addpath('../','-end');  % load main functions (Beamforming)

%% Configuration
maxIter       = 500;  % Number of iterations
modList       = [64 32 16 8 4 2];  % Modulation list
nTxAntennas1  = 1;  % Number of antennas for Radio 1 (192.168.10.1)
nTxAntennas2  = 1;  % Number of antennas for Radio 1 (192.168.110.1)
nTxAntennas   = nTxAntennas1 + nTxAntennas2;
gain          = 8;  % For short distances, use 8. 
NFFT          = 256;  % Do not modify
fileNameReal  = 'weights_real.bin';  % by default, located in /tmp/ dir
fileNameImag  = 'weights_imag.bin';  % by default, located in /tmp/ dir
fileNameTime  = 'weights_time.bin';  % by default, located in /tmp/ dir
flag          = 'real';  % Indicates whether the ch est is real, hardcoded or random
verbose       = false;  % debugging
%% Water-filling configuration
SINRdb        = 20;  % desired SNR
DeltaP        = 0;  % robustness factor (power in linear)
inc           = 0.001;  % water filling increment

%% Configure radios
fprintf('LOG - Configuring experiments...\n');
[radioConfig,bits,~,payload1,trainingSig,fidReal,fidImag,fidTime,channelEst_old] = ...
                f_loadConfig_RX(nTxAntennas,fileNameReal,fileNameImag,fileNameTime);

% Connect to radio
fprintf('LOG - Creating RX System Object. Please wait...\n');
receiver = comm.SDRuReceiver('Platform',radioConfig.rxPlatform1,...
                             radioConfig.rxIDProp1,radioConfig.rxID1,...
                             'MasterClockRate',radioConfig.rxMasterClockRate1,...
                             'DecimationFactor',radioConfig.rxDecimationfactor1,...
                             'ClockSource','External',...
                             'Gain',gain,...
                             'CenterFrequency',900e6,...
                             'EnableBurstMode',true,...
                             'SamplesPerFrame',200e3,...
                             'OutputDataType','double',...
                             'NumFramesInBurst',0);

% Variable to store the BER
BER = zeros(maxIter,length(modList));

%% Main loopif nargin == 2
fprintf('LOG - Begin reception over main loop...\n');
chTot = zeros(nTxAntennas,maxIter);
for i = 1:maxIter
    [rxSig, len] = receiver();
    if len > 0
        %% Frame synchronization + Channel estimation
        [chEst, payload_rx, time_sync] = BER_distributed_chEst(rxSig, trainingSig, [1 1]);

        %% Compute beamforming weights
        BBPowPayload = (payload1'*payload1)/length(payload1);
        w_max = 1/max(real(payload1)) + 1i*1/max(imag(payload1));
        BBPowMax = w_max*w_max';
        [beamWeight,Niter] = f_BFPayload_WF(chEst,DeltaP,inc,SINRdb,0,0,105,BBPowMax,BBPowPayload);

        %% Channel estimation report
        if ~any(isnan(chEst)) && length(chEst)==nTxAntennas
            % Write to file and transmit to TX-hosts using Python - Real
            fseek(fidReal,0,'bof');  % seek begining of the file
            fwrite(fidReal, real(chEst), 'double');
            % Write to file and transmit to TX-hosts using Python - Imag
            fseek(fidImag,0,'bof');  % seek begining of the file
            fwrite(fidImag, imag(chEst), 'double'); 
            % Write to file and transmit to TX-hosts using Python - Time
            fseek(fidTime,0,'bof');  % seek begining of the file
            fwrite(fidTime, time_sync, 'double');
            % Update previous channel estimation
            channelEst_old = chEst;
            % Log channel estimation
            if verbose
                fprintf('Iter %d - ',i);
                for id = 1:nTxAntennas1
                    fprintf('H1(%d) = %.7f + %.7fj | ',id,real(chEst(id)),imag(chEst(id)));
                end
                for id = nTxAntennas1+1:nTxAntennas
                    fprintf('H2(%d) = %.7f + %.7fj | ',id,real(chEst(id)),imag(chEst(id)));
                end
                fprintf('\n');
            end
        end
        chTot(:,i) = channelEst_old;  % Store value in global variable
        
        %% OFDM Demodulaton
        fftOut = fft(reshape(payload_rx, NFFT, 64));

        %% OFDM Subcarrier indexing
%         figure(20); clf('reset');
        for modIdx = 1:length(modList)
            index = 4 + modIdx;
            y = fftOut(index,:).';  % Extract Subcarrier
            y = y/sqrt(mean(y'*y));  % Normalize symbols
            y = 1/sqrt(sum(var(y))).*y;  % Normalize symbols
%             plotMyConstellation(modIdx,modList,bits,y);
            % Compute Bit Error Rate for the 64-QAM modulation
            if ~isempty(y) && ~any(isnan(y))
                % Demodulator expecting normalized symbols
                M = modList(modIdx);
                data_rx = qamdemod(y,M,'OutputType','bit','UnitAveragePower',true);
                % Compute BER
                BER(i,modIdx) = sum(abs(bits{modIdx} - data_rx))/length(data_rx);
                flag = '(real)';
            elseif i>1
                BER(i,modIdx) = BER(i-1,modIdx);  % First element is the BER
                flag = '(Hardcoded)';
            else
                BER(i,modIdx) = 0.5;  % Assigning random value
                flag = '(Random)';
            end
            fprintf('  BER: %.3f - %s (M=%d)\n',BER(i,modIdx),flag,modList(modIdx));
        end
        
    else
        chTot(:,i) = chEst_old;
    end
end

save(strcat('KRIExperiments_',string(datetime)));

lastIter = i;
for idxFFT = 1:nTxAntennas
    figure(15);
    subplot(nTxAntennas,3,3*(idxFFT-1) + 1); hold on; grid minor;
    % plot((1:maxIter),real(chRealTot),'Color','b');
    plot((1:lastIter),real(chTot(idxFFT,1:lastIter)),'Color','r','LineWidth',2);
    title('Real','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
    subplot(nTxAntennas,3,3*(idxFFT-1) + 2); hold on; grid minor;
    % plot((1:maxIter),imag(chRealTot),'Color','b');
    plot((1:lastIter),imag(chTot(idxFFT,1:lastIter)),'Color','r','LineWidth',2);
    title('Imaginary','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);
    subplot(nTxAntennas,3,3*(idxFFT-1) + 3); hold on; grid minor;
    plot((1:lastIter),abs(chTot(idxFFT,1:lastIter)),'Color','r','LineWidth',2);
    title('Absolute Gain','FontSize',12);
    xlabel('Iteration','FontSize',12);
    ylabel('Gain','FontSize',12);

    figure(16);
    subplot(4,1,idxFFT); hold on;
    covTot_imag = zeros(lastIter,1);
    covTot_real = zeros(lastIter,1);
    re = real(chTot(idxFFT,:));
    im = imag(chTot(idxFFT,:));
    re(isnan(re))=0;
    window_size = 1;
    for k = 1:1:lastIter-window_size
        covTot_imag(k) = cov(im(k:k+window_size));
        covTot_real(k) = cov(re(k:k+window_size));
    end
    plot(abs(1-covTot_imag),'LineWidth',2);
    plot(abs(1-covTot_real),'LineWidth',2);
    legend('Imaginary','Real');
    ylim([0.999 1]);
end

figure(17);
totGain = 0;
for txID = 1:nTxAntennas
    t = abs(chTot(txID,1:lastIter));
    t(isnan(t))=0;  % replace nan values
    totGain = totGain + t;
end
 plot((1:lastIter),totGain,'Color','r','LineWidth',2);
title('Total channel Gain','FontSize',12);
xlabel('Iteration','FontSize',12);
ylabel('Gain (linear)','FontSize',12);

figure(18);
subplot(611); plot(abs(BER(1:i,1)),'LineWidth',1.5); legend('64-QAM'); ylabel('BER'); ylim([0 1]);
subplot(612); plot(abs(BER(1:i,2)),'LineWidth',1.5); legend('32-QAM'); ylabel('BER'); ylim([0 1]);
subplot(613); plot(abs(BER(1:i,3)),'LineWidth',1.5); legend('16-QAM'); ylabel('BER'); ylim([0 1]);
subplot(614); plot(abs(BER(1:i,4)),'LineWidth',1.5); legend('8-QAM'); ylabel('BER'); ylim([0 1]);
subplot(615); plot(abs(BER(1:i,5)),'LineWidth',1.5); legend('QPSK'); ylabel('BER'); ylim([0 1]);
subplot(616); plot(abs(BER(1:i,6)),'LineWidth',1.5); legend('BPSK'); ylabel('BER'); ylim([0 1]);




function plotMyConstellation(modIdx,modList,bits,y)
    % Plot constellation
    figure(20);  subplot(1,length(modList),modIdx);  hold on;
    y_tx = qammod(bits{modIdx},modList(modIdx),'InputType','bit','UnitAveragePower',true);
    plot(real(y_tx),imag(y_tx),'LineStyle','None','Marker','.','Color','r');
    plot(real(y),imag(y),'LineStyle','None','Marker','.','Color','b');
    xlim([-2 2]);  ylim([-2 2]);  % Normalized
    tit = strcat('Receiver with k =',{' '},num2str(modList(modIdx)));
    title(tit{1},'FontSize',12);
end