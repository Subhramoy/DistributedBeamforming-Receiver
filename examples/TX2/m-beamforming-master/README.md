# Matlab Beamforming

1. For experiments, only two scripts are relevant for experiments: *BER_centralized_TX* (transmitter) and *BER_centralized_RX* (receiver). The rest are functions called by those.
2. For 1 and 2 antennas, the system uses the TX radio whose IP is *192.168.10.2*. For 3 and 4 antennas, it connects to the additional radio at *192.168.110.2*.

### IMPORTANT: 
The Matlab version should be **R2018a**. It uses the UHD 003.009.007_vendor, which has shown to be more stable than the others. For instance, Matlab R2018b uses a UHD version that fills the buffer and requires a pause in order to update the beamforming weights accordingly. However, this pause generates issues of synchronism between the two daughter-boards in the radio. That is, there's a time gap between the transmissions, which results in ISI at the receiver. Bottom line, use R2018a.

3. The FPGA Firmware at the radios should be updated according to the Matlab R2018a. To do so, just run the matlab function *sdruload* per transmitter (X310). For instance, to update the firmware at the first TX radio we need to call: *sdruload('Device','X310','IPAddress','192.168.10.2')*. The B210 update their firmware automatically during the first Matlab call.
4. When experiencing issues, i.e. high BER, unplug the Octoclock from power source and power cycle the radio that generates the reference signal. Sometimes the Octoclock locks onto an old ref. signal and somehow causes errors in the system.
5. Centralized version means that TX radios and RX should be connected to one computer. Only 2 Matlab instances need to be up. One for TX and the other one for RX.
6. The current beamforming algorithm implemented is a standard channel equalization per-antenna. The future beamforming algorithm should be implemented under the function *f_BFPayload.m*
7. The center frequency should be set at 900MHz. Other frequencies have returned unstable results.