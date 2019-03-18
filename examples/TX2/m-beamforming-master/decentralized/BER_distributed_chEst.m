function [channelEstimate, payload, timeSync] = BER_distributed_chEst(rxSig, goldSeqRef, nTxAntennas)

figure(1); cla reset; hold on;

rxSigLen = length(rxSig);
payloadLen = 64*256;  % 64 symbols with IFFT length of 256

nRadios = length(nTxAntennas);  % determine number of transmitted radios

%% Get gold sequence distribution
goldSeqIdx = getGoldIndices(nTxAntennas);

%% Synchronize and perform channel estimation
crossCorr = cell(nRadios,1);  % store cross-correlation with gold sequences
channelEstimate = zeros(1,goldSeqIdx{end}(end));  % final channel estimation
myStartIdx = cell(nRadios,1);  % store correlation indices
payload = zeros(payloadLen,nRadios);  % store received payload per TX
for RadioID = 1:nRadios
    
    myNTxAntennas = nTxAntennas(RadioID);
    crossCorr{RadioID} = cell(1,myNTxAntennas);
    
    myGoldIdx = goldSeqIdx{RadioID};
    myGoldSeq = goldSeqRef(:,myGoldIdx);
    
    for txID = 1:myNTxAntennas
        %% Detect training signal by cross correlation
        myCorr = xcorr(rxSig,myGoldSeq(:,txID));
        crossCorr{RadioID}{txID} = myCorr((rxSigLen+1):end);
        plot(abs(crossCorr{RadioID}{txID}));
    end
    
    %% Frame synchronization
    [myStartIdx{RadioID},rxGoldSeqRef,payload(:,RadioID)] = ...
            synchronize(crossCorr{RadioID},myGoldSeq,rxSig);

    %% Channel estimation
    chEst = channel_estimation(rxGoldSeqRef,myGoldSeq);
    
    %% Append results
    channelEstimate(1,myGoldIdx) = chEst;
end

%% Compute time syncronization
timeSync = getTimeSync(myStartIdx);

%% Correct received payload
payload = sum(payload,2) ./ nRadios;


% END OF MAIN FUNCTION
end


%% FRAME SYNCHRONIZATION
function [myStartIdx,rxGoldSeqRef,payload] = synchronize(crossCorr,goldSeqRef,rxSig)
    % Compute basics
    nTxAntennas = length(crossCorr);
    rxSigLen = length(rxSig);
    L = size(goldSeqRef,1);
    payloadLen = 64*256;  % 64 symbols with IFFT length of 256
    
    % At least 80% of global maximum
    peakIntervals = find(abs(crossCorr{1})>(0.8*max(abs(crossCorr{1}))));
    
    % We know that the training signal has more than 1000 samples
    peakCandidates1 = find((peakIntervals(2:end)-peakIntervals(1:end-1))>1000);
    if ~isempty(peakCandidates1)
        peakCandidates1(end) = [];  % Don't use the last one
    end
    
    myStartIdx = zeros(1,length(peakCandidates1));
    payload = zeros(payloadLen,1);
    numPayload = 0;
    rxGoldSeqRef = zeros(size(goldSeqRef));
    for i = 1:length(peakCandidates1)
        startIndex = peakIntervals(peakCandidates1(i)+1);
        myStartIdx(i) = startIndex;
        
        % Focus on a small window that contains a peak. +/- 100 samples
        if startIndex > 100; windowRange = startIndex+(-100:100);
        else;                windowRange = 1:(startIndex+100);
        end

        dataWindow = [];
        for txID = 1:nTxAntennas
            dataWindow = [dataWindow crossCorr{txID}(windowRange)];  %#ok<AGROW>
        end
        dataWindow = abs(dataWindow);

        % Find exact locations of the peaks (one in each of 4 cross correlations)
        [~,maxLoc] = max(dataWindow); % maxLoc is a 1x4 vector
        finalPeakLoc = startIndex-101;

        indexAll = finalPeakLoc + maxLoc + rxSigLen;
        indexMin = min(indexAll);
        indexMax = max(indexAll);

        if (indexMax-rxSigLen+L+400+payloadLen) <= rxSigLen % 400 samples between training signal and payload
            % Got a complete payload
            % Extract received training signal and payload
            rxGoldSeqRef = rxSig((indexMin-rxSigLen+1):(indexMax-rxSigLen+L));
            payload = payload + ...
              rxSig((indexMax-rxSigLen+L+400+1):(indexMax-rxSigLen+L+400+payloadLen));
            numPayload = numPayload + 1;
        end
    end
    
    if numPayload > 0
        % The transmitter sends the same payload continously.
        % Return the average of all detected payloads in rxSig.
        payload = payload/numPayload;
    end
end



%% CHANNEL ESTIMATION
function chEst = channel_estimation(rxGoldSeqRef,goldSeqRef)
    % Use least squares fit to estimate channel response
    chEst = goldSeqRef\rxGoldSeqRef;
    chEst = chEst.';
end



%% GOLD SEQUENCES INDICES
function goldSeqIdx = getGoldIndices(nTxAntennas)
    nRadios = length(nTxAntennas);
    CSNAntennas = cumsum(nTxAntennas);  % vector with cum sum of number of antennas
    goldSeqIdx = cell(1,nRadios);
    for RadioID = 1:nRadios
        if RadioID > 1
            goldSeqIdx{RadioID} = ( CSNAntennas(RadioID-1)+1 : CSNAntennas(RadioID) );
        else
            goldSeqIdx{RadioID} = ( 1 : CSNAntennas(RadioID) );
        end
    end
end



%% TIME SYNCHRONIZATION
function timeSync = getTimeSync(myStartIdx)
    nRadios = length(myStartIdx);
    timeSync = zeros(1,nRadios);
    
    if any(cellfun(@isempty,myStartIdx))
        return;
    end
    
    % Compute closest sample difference between radios
    myDiff = zeros(1,nRadios);
    myDiff(1) = 0;  % first radio is sync with itself
    % Get minimum length
    for radioID = 1:nRadios;  myLengths(radioID) = length(myStartIdx{radioID});  end
    [minLen,~] = min(myLengths);
    % Trim indices to have the same length
    for radioID = 1:nRadios;  myStartIdx{radioID} = myStartIdx{radioID}(1:minLen);  end
    
    for radioID = 2:nRadios
        % always take first goldSeq as reference
        substract = myStartIdx{1} - myStartIdx{radioID};
        [~,idx] = min(abs(substract));
        myDiff(radioID) = substract(idx);
    end
    
    % NOTE: A negative difference implies that the first one is ahead
    % NOTE: Time correction applies for those that are ahead
    [~,myTimeRef] = min(myDiff);
    timeSync = myDiff - myDiff(myTimeRef);
end



