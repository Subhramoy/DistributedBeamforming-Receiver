function [beamWeight,Niter] = f_BFPayload_WF(chEst,DeltaP,inc,SINRdb,Atx,Gtx,Noise,BBPowMax,BBPowPayload)

% Initialize parameters
chEst_abs = abs(chEst);  % Channel gain 
w_abs = zeros(1,size(chEst,2));  % Inizialization of weight power
Niter = 1;  % initialize number of total iterations  

% Check if requested SNR is attainable
P_ch = -pow2db(chEst*chEst');
BBPowMax_rep = repmat(BBPowPayload*BBPowMax,1,length(chEst));
P_tx_offered = pow2db(sum(BBPowMax_rep)) + 30;  % in dBm
SNRdb_offered = P_tx_offered + Atx + Gtx - P_ch + Noise;

if SNRdb_offered < SINRdb
    fprintf('WARNING - Cannot achieve SNR. %.3f vs %.3f (dB)\n',SNRdb_offered,SINRdb);
    fprintf('WARNING - Consider revising Gains\n');
    w_abs = sqrt(BBPowMax);  % assign maximum
else
    % Water filling - power allocation
    targetBB = SINRdb + P_ch - Atx - Gtx - Noise - 30;  % in dB (not in dBm, important)
    while pow2db(BBPowPayload*(w_abs*w_abs')) < targetBB
        [~,idxvalids] = find(w_abs.^2 + inc < BBPowMax);
        if isempty(idxvalids)
            break;
        end
        [~, in] = sort(w_abs(idxvalids).^2 + (1./chEst_abs(idxvalids)).^2,'ascend');  % Locate the lowest noise level channel
        w_abs(idxvalids(in(1))) = w_abs(idxvalids(in(1)))+inc;  % Otherwise go to the next smallest power channel
%         fprintf('iter %d - Achieved Power = %.3f (dBm)\n',Niter,pow2db(w_abs*w_abs'));
        Niter = Niter + 1;  % increase number of iterations
    end
    
    % Robust version
    w_abs_max = sqrt(BBPowMax);
    w_abs = min(w_abs + DeltaP,w_abs_max);
end

% Assign weights (phase and power)
beamWeight_angle = (-1.*angle(chEst)).';          % Phase equalization
beamWeight = w_abs.*(cos(beamWeight_angle).'+...  % Complex beamWeights
                       1i.*sin(beamWeight_angle).');



% EOF