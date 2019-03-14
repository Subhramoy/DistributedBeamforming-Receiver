#!/usr/bin/python

import numpy


# Test inputs

channel_est = [1+1j , 2+2j, 3+3j , 4+4j]

def pow2db(a):
    return 10*numpy.log10(a)

def calculate_weights(channel_estimations):
    ## Static variables
    DeltaP = 0
    inc = 0.001
    SINRdb = 20
    Atx = 0
    Gtx = 0
    Noise = 105
    BBPowMax = 6327e+03
    BBPowPayload = 9.1553e-05

    print BBPowPayload
    chEst_abs = numpy.absolute(channel_estimations);  # Channel gain
    w_abs = numpy.zeros(len(channel_estimations), dtype=float);  # Inizialization of weight power
    Niter = 1;  # initialize number of total iterations

    ## Check if requested SNR is attainable
    # P_ch = -pow2db(chEst*chEst');
    P_ch = -1*pow2db( # pow2db
                        numpy.matmul( # Matrix multiply
                            channel_estimations, # CE array
                            numpy.matrix(channel_estimations).getH() # Complex conjugate
                            )
                        )


    #BBPowMax_rep = repmat(BBPowPayload*BBPowMax,1,length(channel_estimations));
    BBPowMax_rep = numpy.repeat(BBPowPayload*BBPowMax,len(channel_estimations));

    #P_tx_offered = pow2db(sum(BBPowMax_rep)) + 30;
    P_tx_offered = pow2db(numpy.sum(BBPowMax_rep)) + 30;

    #SNRdb_offered = P_tx_offered + Atx + Gtx - P_ch + Noise;
    SNRdb_offered = P_tx_offered + Atx + Gtx - P_ch + Noise;

    if SNRdb_offered < SINRdb:
        print('WARNING - Cannot achieve SNR. {}}vs {}} (dB)'.format(SNRdb_offered,SINRdb))
        print('WARNING - Consider revising Gains\n');
        w_abs = sqrt(BBPowMax);  # assign maximum

    else:
        # Water filling - power allocation
        targetBB = SINRdb + P_ch - Atx - Gtx - Noise - 30;  # in dB (not in dBm, important)

        #while pow2db(BBPowPayload*(w_abs*w_abs')) < targetBB
        while pow2db( numpy.matmul(
                                    BBPowPayload,
                                     numpy.matmul(
                                            w_abs,
                                            numpy.matrix(w_abs).getH() # Complex conjugate
                                                )
                                )  < targetBB:

            [~,idxvalids] = find(w_abs.^2 + inc < BBPowMax);
            if isempty(idxvalids)
                break;

            # Locate the lowest noise level channel
            [~, in] = sort(w_abs(idxvalids).^2 + (1./chEst_abs(idxvalids)).^2,'ascend')

            # Otherwise go to the next smallest power channel
            w_abs(idxvalids(in(1))) = w_abs(idxvalids(in(1)))+inc
            print('iter {}- Achieved Power = {} (dBm)\n'.format(Niter,pow2db(w_abs*w_abs')))
            Niter = Niter + 1;

        % Assign weights (phase and power)
beamWeight_angle = (-1.*angle(chEst)).';          % Phase equalization
beamWeight = w_abs.*(cos(beamWeight_angle).'+...  % Complex beamWeights
                       1i.*sin(beamWeight_angle).');



    return BBPowPayload



if __name__ == '__main__':
    print calculate_weights(channel_est)
