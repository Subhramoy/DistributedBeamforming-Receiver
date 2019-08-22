%% Clear environment
clear all; close all; clc;
addpath('utils','-end');

%% Load transmitted signal
load(fullfile('data','information5.mat'),'bits','symbols','payload1');

%% CONFIGURATION
maxIter       = 1000;
fileName      = 'payload.dat';
readSize      = 64*256;  % Payload length is 64*256
modList       = [64 32 16 8 4 2];  % Modulation list
% Create feedback file if inexistent
if ~exist(fullfile(tempdir,fileName),'file')
    % Create feedback file if inexistent
    fid = fopen(fullfile(tempdir,fileName),'wb');
    fclose(fid);
end
% Open same file in reading format
fid = fopen(fullfile(tempdir,fileName),'rb');

%% MAIN LOOP
BER = zeros(maxIter,length(modList));  % Variable to store the BER
for i = 1:maxIter
    %% Check expected length
    myRead = read_complex_binary_CBG(fid, readSize);
    lengthMyRead = length(myRead);
    if lengthMyRead==readSize
        %% OFDM Demodulaton
        fftOut = fft(reshape(myRead, 256, 64));
    
        %% OFDM Subcarrier indexing
         figure(20); clf('reset');
        for modIdx = 1:length(modList)
            index = 4 + modIdx;
            y = fftOut(index,:).';  % Extract Subcarrier
            y = y/sqrt(mean(y'*y));  % Normalize symbols
            y = 1/sqrt(sum(var(y))).*y;  % Normalize symbols
                 plotMyConstellation(modIdx,modList,bits,y);
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
        fprintf('Iter %d - File is empty\n',i);
    end
    pause(0.05);
end

fclose(fid);

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