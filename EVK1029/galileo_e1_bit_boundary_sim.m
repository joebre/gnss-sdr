function galileo_e1_bit_boundary_sim(transition_ms, doppler_step_hz, n_incoherent, f_error_hz, cn0_dBHz)
% GALILEO_E1_BIT_BOUNDARY_SIM  Galileo E1 4 ms coherent + N incoherent
% accumulation across 5 Doppler bins, when the coherent window straddles a
% navigation data bit boundary.
%
% Usage:
%   galileo_e1_bit_boundary_sim()                      % all defaults
%   galileo_e1_bit_boundary_sim(2)                      % bit transition at 2 ms
%   galileo_e1_bit_boundary_sim(2, 125, 5, 40, 35)      % all args
%
% Arguments (all optional, in order; pass [] to keep a default):
%   transition_ms   - time within each 4 ms coherent window (0, 1, 2, or 3)
%                     at which the data bit value changes. 0 means the
%                     window does not straddle a boundary (the whole 4 ms
%                     window is a single bit). Default: 2.
%   doppler_step_hz - spacing between adjacent Doppler bins under test.
%                     Default: 1/(2*T) = 125 Hz for T = 4 ms (matches
%                     Acquisition_1B.doppler_step's usual formula).
%   n_incoherent    - number of incoherent (non-coherent power) dwells
%                     accumulated per bin. Default: 5.
%   f_error_hz      - true residual carrier frequency error of the signal,
%                     relative to the bin grid's center (the "0 Hz" bin).
%                     Default: 40.
%   cn0_dBHz        - carrier-to-noise density ratio, dB-Hz. Default: 35.
%
% Model / simplifications:
%   The PRN code and BOC(1,1)/CBOC subcarrier correlation are assumed
%   perfect -- this isolates the bit-boundary + frequency-error effect
%   only, not a full code-correlation simulation. The signal is an ideal
%   baseband tone at f_error_hz, multiplied by a random +-1 navigation
%   data bit, split into two segments at transition_ms if the window
%   straddles a bit edge. Each of the n_incoherent 4 ms dwells redraws its
%   own random bit pair independently: consecutive real navigation bits
%   are themselves independent random data, so this reproduces the correct
%   per-dwell statistics without needing to chain bit continuity across
%   dwells. Complex AWGN, matching cn0_dBHz, is added at the sample level
%   before correlation.
%
%   For each of the 5 Doppler bins (0, +-1, +-2 times doppler_step_hz),
%   the same n_incoherent dwells are correlated against that bin's
%   frequency hypothesis (unaware of the random bit, exactly as a real
%   acquisition/frequency-scan correlator is) and their power (|.|^2) is
%   summed incoherently. Plots the 5 resulting accumulated powers as a
%   bar chart, and prints the winning bin.

    if nargin < 1 || isempty(transition_ms), transition_ms = 2; end
    if nargin < 2 || isempty(doppler_step_hz), doppler_step_hz = 1 / (2 * 4e-3); end
    if nargin < 3 || isempty(n_incoherent), n_incoherent = 5; end
    if nargin < 4 || isempty(f_error_hz), f_error_hz = 40; end
    if nargin < 5 || isempty(cn0_dBHz), cn0_dBHz = 35; end

    if ~ismember(transition_ms, [0 1 2 3])
        error('transition_ms must be one of 0, 1, 2, 3');
    end

    T = 4e-3;               % coherent integration time [s]
    fs = 20000;             % simulation sample rate [Hz] (baseband model only)
    A = 1;                  % normalized signal amplitude

    % C/N0 -> per-sample complex AWGN std. N0 = noise_power / fs, so
    % C/N0 = A^2 * fs / noise_power  =>  noise_power = A^2 * fs / (C/N0).
    cn0_linear = 10 ^ (cn0_dBHz / 10);
    noise_power = (A ^ 2) * fs / cn0_linear;
    noise_sigma = sqrt(noise_power);   % total complex noise std (I and Q each noise_sigma/sqrt(2))

    bin_multipliers = [-2 -1 0 1 2];
    bin_freqs_hz = bin_multipliers * doppler_step_hz;
    n_bins = numel(bin_freqs_hz);

    accumulated_power = zeros(1, n_bins);

    t1 = transition_ms * 1e-3;   % duration of the leading bit segment
    t2 = T - t1;                 % duration of the trailing bit segment
    n1 = round(t1 * fs);
    n2 = round(t2 * fs);

    for dwell = 1:n_incoherent
        % Random +-1 navigation data bits for this dwell's two segments.
        b1 = 2 * (rand() > 0.5) - 1;
        b2 = 2 * (rand() > 0.5) - 1;

        if n1 > 0
            t_seg1 = (0:n1 - 1) / fs;
            sig1 = A * b1 * exp(1j * 2 * pi * f_error_hz * t_seg1);
        else
            sig1 = [];
        end
        if n2 > 0
            t_seg2 = (0:n2 - 1) / fs + t1;
            sig2 = A * b2 * exp(1j * 2 * pi * f_error_hz * t_seg2);
        else
            sig2 = [];
        end
        sig = [sig1, sig2];
        n_samples = numel(sig);

        % Same noisy input samples are tested against every bin this dwell,
        % exactly as a real Doppler-bin scan reuses one set of input samples.
        noise = (noise_sigma / sqrt(2)) * (randn(1, n_samples) + 1j * randn(1, n_samples));
        rx = sig + noise;
        t_full = (0:n_samples - 1) / fs;

        for k = 1:n_bins
            local_replica = exp(-1j * 2 * pi * bin_freqs_hz(k) * t_full);
            corr = sum(rx .* local_replica);
            accumulated_power(k) = accumulated_power(k) + abs(corr) ^ 2;
        end
    end

    figure;
    bar(bin_freqs_hz, accumulated_power);
    xlabel('Doppler bin offset from grid center [Hz]');
    ylabel('Accumulated power (|.|^2, summed over incoherent dwells)');
    title(sprintf(['Galileo E1 bit-boundary sim: transition=%d ms, step=%.1f Hz, ' ...
                   'N_{incoh}=%d, f_{error}=%.1f Hz, CN0=%.1f dB-Hz'], ...
                  transition_ms, doppler_step_hz, n_incoherent, f_error_hz, cn0_dBHz));
    grid on;
    set(gca, 'XTick', bin_freqs_hz);

    fprintf('\nBin (Hz)      Power\n');
    for k = 1:n_bins
        fprintf('%8.1f   %10.4f\n', bin_freqs_hz(k), accumulated_power(k));
    end
    [~, best_idx] = max(accumulated_power);
    fprintf('\nWinning bin: %.1f Hz (true residual error was %.1f Hz)\n', ...
            bin_freqs_hz(best_idx), f_error_hz);
end
