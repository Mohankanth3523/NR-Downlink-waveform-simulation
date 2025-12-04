%% Multi-Numerology OFDM Simulation (No 5G Toolbox)
% Numerologies: 15, 30, 60, 120 kHz
% Figures:
%  1: 15 kHz spectrum
%  2: 30 kHz spectrum
%  3: 60 kHz spectrum
%  4: 120 kHz spectrum
%  5: Combined spectrum (all SCS)
%  6: Multipath channel (power delay profile)
%  7: BER vs SNR
%  8: EVM vs SNR
%  9: Throughput vs SNR

clear; clc; close all;

%% Simulation parameters
SCS_list_kHz   = [15 30 60 120];   % 4 numerologies
SNR_dB_vec     = 0:5:30;           % SNR points

Nfft           = 1024;             % FFT size
Nused          = 600;              % active subcarriers (centered)
Nsym_per_slot  = 14;               % OFDM symbols per slot
numSlots       = 5;                % number of slots
modOrder       = 4;                % QPSK
bitsPerSym     = log2(modOrder);
cpFrac         = 0.07;             % CP fraction

% Simple multipath channel profile
pathDelays_us  = [0 0.3 0.9 2.5];  % in microseconds
pathPowers_dB  = [0 -2 -6 -10];    % relative powers

numSCS = numel(SCS_list_kHz);
numSNR = numel(SNR_dB_vec);

% Result arrays
BER         = zeros(numSCS, numSNR);
EVM_rms_pct = zeros(numSCS, numSNR);
Throughput  = zeros(numSCS, numSNR);

% For spectrum plots
specFreq  = cell(numSCS,1);
specMagdB = cell(numSCS,1);

% For multipath channel visualization (store first numerology's channel)
h_first   = [];
delay_first_us = [];

%% QPSK helpers
qpsk_map = @(bits2) ((1 - 2*bits2(:,1)) + 1j*(1 - 2*bits2(:,2)))/sqrt(2);
qpsk_demap = @(sym) [real(sym) < 0, imag(sym) < 0];

%% Main loop over numerologies
for iSCS = 1:numSCS
    
    scs_kHz = SCS_list_kHz(iSCS);
    fprintf('\n=== SCS = %d kHz ===\n', scs_kHz);
    
    % Subcarrier spacing and OFDM timing
    deltaF  = scs_kHz * 1e3;       % Hz
    Fs      = Nfft * deltaF;       % sampling rate
    Tu      = 1 / deltaF;          % useful symbol duration
    Tcp     = cpFrac * Tu;         % CP duration
    Ncp     = round(Tcp * Fs);     % CP samples
    Ns      = Nfft + Ncp;          % samples per OFDM symbol
    
    Tsymbol = Tu + Tcp;
    Tslot   = Nsym_per_slot * Tsymbol;
    frameDuration = numSlots * Tslot;  % total simulated time
    
    % Active subcarriers (centered)
    k0 = Nfft/2 - Nused/2 + 1;
    k1 = Nfft/2 + Nused/2;
    dataIdx = k0:k1;
    
    % Total OFDM symbols and bits
    numSymbolsTotal = Nsym_per_slot * numSlots;
    numBitsTotal    = Nused * numSymbolsTotal * bitsPerSym;
    
    %% Build discrete multipath channel for this Fs
    pathPowersLin = 10.^(pathPowers_dB/10);
    pathPowersLin = pathPowersLin / sum(pathPowersLin);
    
    pathDelays_s  = pathDelays_us * 1e-6;
    pathDelays_smpl = round(pathDelays_s * Fs);
    
    Lh = max(pathDelays_smpl) + 1;
    h  = zeros(Lh,1);
    for p = 1:length(pathDelays_smpl)
        h(pathDelays_smpl(p)+1) = sqrt(pathPowersLin(p)) .* ...
                                  (randn + 1j*randn)/sqrt(2);
    end
    
    % Store first numerology's channel for multipath figure
    if iSCS == 1
        h_first = h;
        delay_first_us = (0:Lh-1)/Fs * 1e6;   % convert sample index to microseconds
    end
    
    % Frequency response for equalization
    Hf = fftshift(fft(h, Nfft));
    
    %% --- One clean TX waveform for this numerology (for spectrum only) ---
    txBits_spec = randi([0 1], numBitsTotal, 1);
    txBits2_spec = reshape(txBits_spec, 2, []).';
    txSymbols_all_spec = qpsk_map(txBits2_spec);
    txGridData_spec = reshape(txSymbols_all_spec, Nused, numSymbolsTotal);
    
    txWaveform = zeros(numSymbolsTotal*Ns,1);
    for nSym = 1:numSymbolsTotal
        X = zeros(Nfft,1);
        X(dataIdx) = txGridData_spec(:,nSym);
        x    = ifft(ifftshift(X));          % time-domain symbol
        x_cp = [x(end-Ncp+1:end); x];       % add CP
        idx_start = (nSym-1)*Ns + 1;
        idx_end   = nSym*Ns;
        txWaveform(idx_start:idx_end) = x_cp;
    end
    
    % Spectrum calculation
    Nspec = 16384;
    if length(txWaveform) < Nspec
        xSeg = txWaveform;
        Nspec = length(txWaveform);
    else
        xSeg = txWaveform(1:Nspec);
    end
    Xf = fftshift(fft(xSeg, Nspec));
    psd = 20*log10(abs(Xf)/max(abs(Xf) + 1e-12));
    fAxis = (-Nspec/2:Nspec/2-1) * (Fs/Nspec) / 1e6;  % MHz
    
    specFreq{iSCS}  = fAxis;
    specMagdB{iSCS} = psd;
    
    % ---------- Individual spectrum figure for each SCS ----------
    figure(iSCS);  % Figure 1,2,3,4 correspond to 15,30,60,120 kHz
    plot(fAxis, psd, 'LineWidth', 1.5);
    grid on;
    xlabel('Frequency (MHz)');
    ylabel('Normalized Spectrum (dB)');
    title(sprintf('Transmit OFDM Spectrum - %d kHz SCS', scs_kHz));
    
    %% --- Now SNR loop for BER / EVM / Throughput ---
    for iSNR = 1:numSNR
        
        SNRdB = SNR_dB_vec(iSNR);
        fprintf('  SNR = %d dB\n', SNRdB);
        
        % Transmitter: bits -> QPSK -> OFDM
        txBits = randi([0 1], numBitsTotal, 1);
        txBits2 = reshape(txBits, 2, []).';
        txSymbols_all = qpsk_map(txBits2);
        txGridData = reshape(txSymbols_all, Nused, numSymbolsTotal);
        
        txWaveform = zeros(numSymbolsTotal*Ns,1);
        for nSym = 1:numSymbolsTotal
            X = zeros(Nfft,1);
            X(dataIdx) = txGridData(:,nSym);
            x    = ifft(ifftshift(X));
            x_cp = [x(end-Ncp+1:end); x];
            idx_start = (nSym-1)*Ns + 1;
            idx_end   = nSym*Ns;
            txWaveform(idx_start:idx_end) = x_cp;
        end
        
        % Channel + AWGN
        rxChan = conv(txWaveform, h);
        sigPower   = mean(abs(rxChan).^2);
        noisePower = sigPower / (10^(SNRdB/10));
        noise = sqrt(noisePower/2) * (randn(size(rxChan)) + 1j*randn(size(rxChan)));
        rxNoisy = rxChan + noise;
        
        % Align received signal (remove initial channel transient)
        rxAligned = rxNoisy(Lh:end);
        totalSamplesNeeded = numSymbolsTotal * Ns;
        if length(rxAligned) < totalSamplesNeeded
            rxAligned(end+1:totalSamplesNeeded) = 0;
        else
            rxAligned = rxAligned(1:totalSamplesNeeded);
        end
        
        % Receiver: remove CP, FFT, equalize
        rxGridEq = zeros(Nused, numSymbolsTotal);
        for nSym = 1:numSymbolsTotal
            idx_start = (nSym-1)*Ns + 1;
            idx_end   = nSym*Ns;
            y_cp = rxAligned(idx_start:idx_end);
            y    = y_cp(Ncp+1:end);         % remove CP
            Y    = fftshift(fft(y, Nfft));  % frequency-domain
            Yd   = Y(dataIdx);
            Hd   = Hf(dataIdx);
            eqSymbols = Yd ./ (Hd + 1e-12); % 1-tap equalizer
            rxGridEq(:, nSym) = eqSymbols;
        end
        
        rxSymbols_all = rxGridEq(:);
        
        % --------- EVM ----------
        Nmin  = min(length(txSymbols_all), length(rxSymbols_all));
        refSym = txSymbols_all(1:Nmin);
        rxSym  = rxSymbols_all(1:Nmin);
        errSym = rxSym - refSym;
        EVM_rms = sqrt(mean(abs(errSym).^2)) / sqrt(mean(abs(refSym).^2));
        EVM_rms_pct(iSCS, iSNR) = 100 * EVM_rms;
        
        % --------- BER ----------
        rxBits2_hat = qpsk_demap(rxSym);
        rxBits_hat  = reshape(rxBits2_hat.', [], 1);
        Nbits_min   = min(length(txBits), length(rxBits_hat));
        txBits_use  = txBits(1:Nbits_min);
        rxBits_use  = rxBits_hat(1:Nbits_min);
        numErrBits  = sum(txBits_use ~= rxBits_use);
        BER(iSCS, iSNR) = numErrBits / Nbits_min;
        
        % --------- Throughput (bits/s) ----------
        numCorrectBits = Nbits_min - numErrBits;
        Throughput(iSCS, iSNR) = numCorrectBits / frameDuration;
        
    end
end

%% 5) Combined spectrum figure (all numerologies)
figure(5);
for iSCS = 1:numSCS
    plot(specFreq{iSCS}, specMagdB{iSCS}, 'LineWidth', 1.5); hold on;
end
grid on;
xlabel('Frequency (MHz)');
ylabel('Normalized Spectrum (dB)');
legend(arrayfun(@(x) sprintf('%d kHz SCS', x), SCS_list_kHz, 'UniformOutput', false), ...
       'Location','southwest');
title('Transmit OFDM Spectrum - All Subcarrier Spacings');

%% 6) Multipath environment figure (power delay profile)
figure(6);
stem(delay_first_us, abs(h_first).^2, 'filled');
grid on;
xlabel('Delay (\mus)');
ylabel('Power |h(\tau)|^2');
title('Simulated Multipath Channel - Power Delay Profile');

%% 7) BER vs SNR
figure(7);
for iSCS = 1:numSCS
    semilogy(SNR_dB_vec, BER(iSCS,:), '-o', 'LineWidth', 1.5); hold on;
end
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend(arrayfun(@(x) sprintf('%d kHz', x), SCS_list_kHz, 'UniformOutput', false), ...
       'Location','southwest');
title('BER vs SNR for Different Subcarrier Spacings');

%% 8) EVM vs SNR
figure(8);
for iSCS = 1:numSCS
    plot(SNR_dB_vec, EVM_rms_pct(iSCS,:), '-o', 'LineWidth', 1.5); hold on;
end
grid on;
xlabel('SNR (dB)');
ylabel('EVM_{RMS} (%)');
legend(arrayfun(@(x) sprintf('%d kHz', x), SCS_list_kHz, 'UniformOutput', false), ...
       'Location','northeast');
title('EVM vs SNR for Different Subcarrier Spacings');

%% 9) Throughput vs SNR
figure(9);
for iSCS = 1:numSCS
    plot(SNR_dB_vec, Throughput(iSCS,:)/1e6, '-o', 'LineWidth', 1.5); hold on;
end
grid on;
xlabel('SNR (dB)');
ylabel('Throughput (Mbit/s)');
legend(arrayfun(@(x) sprintf('%d kHz', x), SCS_list_kHz, 'UniformOutput', false), ...
       'Location','southeast');
title('Throughput vs SNR for Different Subcarrier Spacings');

disp('Simulation complete with separate and combined spectrum + multipath figure.');
