function beamWeight = f_BFPayload(channelEst)

% % Rotate payload 1 so that receiver 1 does not need to
% % correct the phase of payload 1
% beamWeight = ones(size(channelEst));  % Does not matter
% phaseCorrection1 = beamWeight  * channelEst.';  % Beamforming info
% phaseCorrection1 = phaseCorrection1/abs(phaseCorrection1);
% beamWeight = beamWeight ./ phaseCorrection1;

beamWeight = 1./channelEst;
beamWeight = beamWeight./max(abs(beamWeight));

% EOF