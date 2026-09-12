% fll_discriminator_sim.m
%
% Open-loop (no tracking loop) comparison of two FLL discriminators, both
% evaluated on noisy correlator samples for a given carrier-to-noise density
% ratio (CN0) and a given true carrier frequency error:
%
%   1) "Discriminator-P" (new): the energy-based FLL discriminator from
%      Xinhua Tang's PhD thesis, "Development and Analysis of Advanced
%      Techniques for GNSS Receivers" (Politecnico di Torino, 2014), Ch. 4,
%      "New Design of a FLL Discriminator Based on Energy". It needs three
%      Prompt correlations per epoch, taken with local carrier replicas at
%      f_P-f_step, f_P, f_P+f_step (f_step = 2/(3*Td)):
%          f_error = (|Z_R|-|Z_L|) / (|Z_L|+|Z_P|+|Z_R|) * f_step
%      (Sign flipped from the thesis' literal Eq. (4.9) -- see NOTE below.)
%
%   2) The traditional four-quadrant ATAN FLL discriminator, using the
%      Prompt correlation from two consecutive epochs:
%          f_error = atan2(cross,dot) / (2*pi*Td)
%
% This mirrors the C++ discriminators added to GNSS-SDR's
% dll_pll_veml_tracking.cc / tracking_discriminators.cc for the EVK1029
% Galileo E1 tracking experiment (Tracking_1B.enable_energy_fdiscriminator).
% No closed tracking loop is simulated here -- only the discriminators'
% raw input/output behavior under noise, i.e. their (noisy) S-curve.
%
% NOTE on sign convention: the thesis' literal Eq. (4.9) reads
%     f_error = (|Z_L|-|Z_R|) / (|Z_L|+|Z_P|+|Z_R|) * f_step
% Under the Z_L=f_P-f_step / Z_R=f_P+f_step correlator assignment (Fig. 4.1),
% that expression's S-curve has NEGATIVE slope through zero, i.e. it points
% the wrong way for closed-loop negative feedback. Confirmed analytically
% (raising the local frequency reduces the true-vs-local mismatch and
% increases |Z|, so |Z_R|>|Z_L| for a positive true-minus-local error) and
% empirically (the literal formula caused far more frequent loss-of-lock
% when driving GNSS-SDR's tracking loop against a real EVK1029 Galileo E1
% capture). The sign-corrected form used here matches the increasing S-curve
% shown in the thesis' own Fig. 4.2.
%
% Usage:
%   fll_discriminator_sim()                 % defaults: CN0=35 dB-Hz, f_error=50 Hz
%   fll_discriminator_sim(35, 50)           % CN0, true frequency error
%   fll_discriminator_sim(25, 200, 4e-3, 5000, 1)   % + Td, n_trials, seed
%
% NOTE on the default frequency error: the energy discriminator's Left/Right
% branches have a sinc-shaped envelope that goes through zero (a "null")
% right at the edge of its nominal linear range, +/-f_step/2 (e.g. +/-83.3 Hz
% for Td=4 ms). Near that null the branch's noiseless signal component
% nearly vanishes, so its |Z| becomes noise-dominated (Rician-with-small-mean
% behavior) and the discriminator output develops extra bias/variance -- a
% real, expected property of this discriminator, not a simulation bug. The
% default f_error_hz=50 Hz sits comfortably inside the linear range to
% demonstrate typical operation; try values near/beyond f_step/2 to see that
% degradation for yourself.
function fll_discriminator_sim(CN0_dBHz, f_error_hz, Td, n_trials, seed)

    %% ------------------------- PARAMETERS --------------------------------
    if nargin < 1 || isempty(CN0_dBHz),  CN0_dBHz  = 30;    end  % carrier-to-noise density ratio [dB-Hz]
    if nargin < 2 || isempty(f_error_hz), f_error_hz = 50;  end  % true carrier frequency error [Hz]
    if nargin < 3 || isempty(Td),        Td        = 4e-3;  end  % coherent integration time [s]
    if nargin < 4 || isempty(n_trials),  n_trials  = 5000;  end  % Monte Carlo trials at the operating point
    if nargin < 5 || isempty(seed),      seed      = 1;     end  % RNG seed, [] = random
    %% ----------------------------------------------------------------------

    if ~isempty(seed)
        rand('state', seed);
        randn('state', seed);
    end

    f_step = 2 / (3 * Td);                       % thesis Eq. (4.8)
    A      = sqrt(2 * 10^(CN0_dBHz / 10) * Td);  % correlator amplitude: SNR = 2*(C/N0)*Td, unit-variance complex noise per branch

    %% -------- Monte Carlo at the single operating point f_error_hz --------
    out_energy = zeros(1, n_trials);
    out_atan   = zeros(1, n_trials);
    for i = 1:n_trials
        [ZL, ZP, ZR]  = gen_correlations(A, f_error_hz, Td, f_step);
        [Zprev, Zcur] = gen_atan_pair(A, f_error_hz, Td);

        out_energy(i) = disc_energy(ZL, ZP, ZR, f_step);
        out_atan(i)   = disc_atan(Zprev, Zcur, Td);
    end

    printf('--- fll_discriminator_sim: CN0=%.1f dB-Hz, true f_error=%.1f Hz, Td=%.1f ms, %d trials ---\n', ...
           CN0_dBHz, f_error_hz, Td * 1e3, n_trials);
    printf('                    mean [Hz]   std [Hz]   RMSE vs truth [Hz]\n');
    printf('  Energy disc.      %9.3f  %9.3f  %9.3f\n', ...
           mean(out_energy), std(out_energy), sqrt(mean((out_energy - f_error_hz).^2)));
    printf('  Traditional ATAN  %9.3f  %9.3f  %9.3f\n', ...
           mean(out_atan), std(out_atan), sqrt(mean((out_atan - f_error_hz).^2)));

    figure('Name', 'FLL discriminator output distribution');
    nbins = 60;
    subplot(2, 1, 1);
    hist(out_energy, nbins);
    hold on;
    yl = ylim();
    line([f_error_hz, f_error_hz], yl, 'Color', 'r', 'LineWidth', 1.5);
    xlabel('Discriminator output [Hz]'); ylabel('Count');
    title(sprintf('Energy discriminator (CN0=%.0f dB-Hz, true error=%.0f Hz)', CN0_dBHz, f_error_hz));
    legend('samples', 'true error', 'Location', 'northeast');

    subplot(2, 1, 2);
    hist(out_atan, nbins);
    hold on;
    yl = ylim();
    line([f_error_hz, f_error_hz], yl, 'Color', 'r', 'LineWidth', 1.5);
    xlabel('Discriminator output [Hz]'); ylabel('Count');
    title('Traditional ATAN discriminator');
    legend('samples', 'true error', 'Location', 'northeast');

    %% ------------------- Noisy S-curve sweep, same CN0 ---------------------
    f_sweep     = linspace(-2 * f_step, 2 * f_step, 41);
    n_sweep_trials = max(200, round(n_trials / 10));
    mean_energy = zeros(size(f_sweep));
    std_energy  = zeros(size(f_sweep));
    mean_atan   = zeros(size(f_sweep));
    std_atan    = zeros(size(f_sweep));

    for k = 1:numel(f_sweep)
        e = f_sweep(k);
        oe = zeros(1, n_sweep_trials);
        oa = zeros(1, n_sweep_trials);
        for i = 1:n_sweep_trials
            [ZL, ZP, ZR]  = gen_correlations(A, e, Td, f_step);
            [Zprev, Zcur] = gen_atan_pair(A, e, Td);
            oe(i) = disc_energy(ZL, ZP, ZR, f_step);
            oa(i) = disc_atan(Zprev, Zcur, Td);
        end
        mean_energy(k) = mean(oe); std_energy(k) = std(oe);
        mean_atan(k)   = mean(oa); std_atan(k)   = std(oa);
    end

    figure('Name', 'FLL discriminator S-curve under noise');
    plot(f_sweep, f_sweep, 'k:', 'LineWidth', 1); hold on;
    he = errorbar(f_sweep, mean_energy, std_energy, 'b-');
    ha = errorbar(f_sweep, mean_atan, std_atan, 'g-');
    set(he, 'LineWidth', 1.3);
    set(ha, 'LineWidth', 1.3);
    grid on;
    xlabel('True input frequency error [Hz]');
    ylabel('Discriminator output [Hz] (mean +/- std)');
    title(sprintf('Discriminator output vs. true error (CN0=%.0f dB-Hz, Td=%.1f ms, %d trials/point)', ...
                  CN0_dBHz, Td * 1e3, n_sweep_trials));
    legend('ideal (y=x)', 'Energy discriminator', 'Traditional ATAN', 'Location', 'northwest');
    yl = ylim();
    line([-f_step / 2, -f_step / 2], yl, 'Color', [0.6 0.6 0.6], 'LineStyle', ':');
    line([f_step / 2, f_step / 2], yl, 'Color', [0.6 0.6 0.6], 'LineStyle', ':');
end


%% ============================= Helper functions ==========================

function y = sinc0(x)
    % Unnormalized sinc: sin(x)/x, with sinc0(0) = 1.
    if x == 0
        y = 1;
    else
        y = sin(x) / x;
    end
end


function [ZL, ZP, ZR] = gen_correlations(A, e_hz, Td, f_step)
    % One noisy epoch's Left/Prompt/Right correlator outputs (thesis Eq.
    % (4.4)), for a residual frequency error e_hz = true - local NCO [Hz].
    % Complex noise has unit variance per I/Q component per branch. The
    % absolute carrier phase is arbitrary (taken as zero) since the energy
    % discriminator only uses |Z|.
    thetaP = e_hz * pi * Td;
    thetaL = (e_hz + f_step) * pi * Td;   % local carrier lowered by f_step -> mismatch increases
    thetaR = (e_hz - f_step) * pi * Td;   % local carrier raised by f_step  -> mismatch decreases

    ZP = A * sinc0(thetaP) * exp(1j * thetaP) + (randn + 1j * randn);
    ZL = A * sinc0(thetaL) * exp(1j * thetaL) + (randn + 1j * randn);
    ZR = A * sinc0(thetaR) * exp(1j * thetaR) + (randn + 1j * randn);
end


function [Zprev, Zcur] = gen_atan_pair(A, e_hz, Td)
    % Noisy Prompt correlations from two consecutive epochs at a constant
    % residual frequency error e_hz, i.e. Zcur's phase has advanced by
    % 2*pi*e_hz*Td relative to Zprev -- the minimum input the traditional
    % ATAN FLL discriminator needs to produce a frequency estimate.
    theta = e_hz * pi * Td;
    Zprev = A * sinc0(theta) * exp(1j * theta)       + (randn + 1j * randn);
    Zcur  = A * sinc0(theta) * exp(1j * (3 * theta)) + (randn + 1j * randn);
end


function fe_hz = disc_energy(ZL, ZP, ZR, f_step)
    % Energy-based FLL discriminator (sign-corrected, see file header note).
    % Independent of signal amplitude and code delay error.
    denom = abs(ZL) + abs(ZP) + abs(ZR);
    if denom == 0
        fe_hz = 0;
        return;
    end
    fe_hz = (abs(ZR) - abs(ZL)) / denom * f_step;
end


function fe_hz = disc_atan(Zprev, Zcurr, Td)
    % Traditional four-quadrant ATAN FLL discriminator (atan2(cross,dot)).
    dotp   = real(Zprev) * real(Zcurr) + imag(Zprev) * imag(Zcurr);
    crossp = real(Zprev) * imag(Zcurr) - real(Zcurr) * imag(Zprev);
    fe_hz  = atan2(crossp, dotp) / (2 * pi * Td);
end
