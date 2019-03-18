function [radioConfig,bits,symbols,payload1,trainingSig,fidReal,fidImag,fidTime,chEst_prelim] = ...
                    f_loadConfig_RX(numTxAntennas,fileNameReal,fileNameImag,fileNameTime)

%% RADIO ENVIRONMENT
if exist(fullfile('data','radioConfig.mat'),'file')
    % Load existing Radio configuration if already configured
    load(fullfile('data','radioConfig.mat'),'radioConfig');
else
    % Find attached radios and allocate radios for transmitter and receivers.
    radioConfig = helperMUBeamformAllocateRadios;
    % Save radio configuration to a file for future use (and receivers)
    save(fullfile('data','radioConfig.mat'),'radioConfig');
end

%% RADIO CONFIGURATION 
% Configure transmitters
radioConfig.rxID1 = '316E259';  % Use selected B210

%% GOLD SEQUENCES (Detection and channel estimation)
if exist(fullfile('data','trainingSig.mat'),'file')
    % Load pre-defined Gold sequences
    load(fullfile('data','trainingSig.mat'),'trainingSig');
else
    % Generate Gold sequences to distribute across users
    trainingSig = helperMUBeamformInitGoldSeq;
    % Save Gold sequences to be used in the future
    save(fullfile('data','trainingSig.mat'),'trainingSig');
end
trainingSig = trainingSig(:,1:numTxAntennas);  % Store only the ones being used

%% PAYLOAD
if exist(fullfile('data','information4.mat'),'file')
    load(fullfile('data','information4.mat'),'bits','symbols','payload1');
else
    % Create transmitted bits over the air
    M = 64;  k = log2(M);  % 64QAM
    bits{1} = randi([0 1],64*k,1);
    M = 32;  k = log2(M);  % 32QAM
    bits{2} = randi([0 1],64*k,1);
    M = 16;  k = log2(M);  % 16QAM
    bits{3} = randi([0 1],64*k,1);
    M = 8;  k = log2(M);  % 8QAM
    bits{4} = randi([0 1],64*k,1);
    M = 4;  k = log2(M);  % QPSK
    bits{5} = randi([0 1],64*k,1);
    M = 2;  k = log2(M);  % BPSK
    bits{6} = randi([0 1],64*k,1);
    M = 4;  k = log2(M);  % QPSK extra
    bits_QPSK_extra = randi([0 1],64*k,2);
    
    % Create modulated symbols
    symbols{1} = qammod(bits{1},64,'InputType','bit','UnitAveragePower',true).';  % 64QAM
    symbols{2} = qammod(bits{2},32,'InputType','bit','UnitAveragePower',true).';  % 32QAM
    symbols{3} = qammod(bits{3},16,'InputType','bit','UnitAveragePower',true).';  % 16QAM
    symbols{4} = qammod(bits{4},8,'InputType','bit','UnitAveragePower',true).';  % 8QAM
    symbols{5} = qammod(bits{5},4,'InputType','bit','UnitAveragePower',true).';  % QPSK
    symbols{6} = qammod(bits{6},2,'InputType','bit','UnitAveragePower',true).';  % BPSK
    sym_QPSK_extra = qammod(bits_QPSK_extra,4,'InputType','bit','UnitAveragePower',true).';  % QPSK extra
    
    % Construct payload 1:
    % 64 symbols with IFFT length of 256. Each symbol uses 8 subcarriers
    % Subcarrier 5 uses 64-QAM. Each point of 64-QAM is used once. Other
    % subcarriers use QPSK
    modOut1 = [zeros(4,64);       % Zero padding
               symbols{1};        % 64-QAM - 1 subcarrier
               symbols{2};        % 32-QAM - 1 subcarrier
               symbols{3};        % 16-QAM - 1 subcarrier
               symbols{4};        % 8-QAM  - 1 subcarrier
               symbols{5};        % QPSK   - 1 subcarrier
               symbols{6};        % BPSK   - 1 subcarrier
               sym_QPSK_extra;    % QPSK   - 1 subcarrier - Extra
               zeros(256-4-8,64)]; % Zero padding
    payload = reshape(ifft(modOut1),[],1);
    % Scale time-domain signal appropriately
    payload = payload/max(real(payload))*0.5;
    
    % Save information
    save(fullfile('data','information.mat'),'bits','symbols','payload');
end

%% FEEDBACK
% Generate initial channel estimation (in case empty file)
chEst_prelim = rand(1,numTxAntennas) + 1j*rand(1,numTxAntennas);
% Open channel estimate file for RX
fidReal = fopen(fullfile(tempdir,fileNameReal),'wb');
fidImag = fopen(fullfile(tempdir,fileNameImag),'wb');
% Open time synchronization file for RX
fidTime = fopen(fullfile(tempdir,fileNameTime),'wb');
% Initialize beamforming weights for RX
fwrite(fidReal,real(chEst_prelim),'double');
fwrite(fidImag,real(chEst_prelim),'double');




% EOF