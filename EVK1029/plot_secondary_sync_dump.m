% plot_secondary_sync_dump.m
%
% Plots one secondary-code sync stall's corr_value time series (as loaded by
% load_secondary_sync_dump), with the +-lock threshold marked, so you can see
% whether correlation trends toward lock, oscillates near zero, or stays flat --
% rather than just its aggregate distribution.
%
% Usage:
%   stalls = load_secondary_sync_dump();
%   plot_secondary_sync_dump(stalls, 1);        % plot stalls(1)
%
% See also: load_secondary_sync_dump
%
function plot_secondary_sync_dump(stalls, idx)
    if idx < 1 || idx > numel(stalls)
        error('plot_secondary_sync_dump: idx %d out of range (1..%d)', idx, numel(stalls));
    end
    s = stalls(idx);

    figure();
    plot(1:numel(s.corr_value), s.corr_value, '-');
    hold on;
    L = s.secondary_code_length;
    plot([1, numel(s.corr_value)], [L, L], 'r--');
    plot([1, numel(s.corr_value)], [-L, -L], 'r--');
    ylim([-L - 2, L + 2]);
    xlabel('acquire\_secondary() call index');
    ylabel('corr\_value');
    title(sprintf('PRN%d ch%d scan %d: corr\\_value over %d calls (lock at +-%d)', ...
        s.prn, s.channel, s.scan_id, numel(s.corr_value), L));
end
