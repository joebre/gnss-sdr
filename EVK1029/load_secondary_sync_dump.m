% load_secondary_sync_dump.m
%
% Loads the CSV written by dll_pll_veml_tracking's secondary-code sync diagnostic
% (dump_secondary_sync_csv(), fired when the bit_synchronization_time_limit_s
% watchdog kills a channel) into an Octave struct array, one element per stalled
% pull-in attempt. Each row of the CSV is one acquire_secondary() call, in the
% order it happened; this groups rows by scan_id.
%
% CSV columns (header row is skipped):
%   scan_id, prn, system_char, channel, secondary_code_length, call_index, corr_value
%
% corr_value is a +-1 running tally over the current 25-symbol (for Galileo E1C)
% sliding window: +1 per symbol matching the known secondary code polarity, -1 per
% mismatch. Lock requires corr_value = +-secondary_code_length exactly (zero
% tolerance for a wrong symbol anywhere in the window).
%
% Usage:
%   stalls = load_secondary_sync_dump();                    % default: ./secondary_sync_dump.csv
%   stalls = load_secondary_sync_dump('/path/to/dump.csv');
%   load_secondary_sync_dump();                              % no output arg: also prints a summary table
%
% Each stalls(i) has fields:
%   scan_id, prn, system_char, channel, secondary_code_length,
%   corr_value (row vector, in call order -- plot this directly, e.g. plot(stalls(i).corr_value))
%
% See also: plot_secondary_sync_dump
%
function stalls = load_secondary_sync_dump(filename)
    if nargin < 1 || isempty(filename)
        filename = 'secondary_sync_dump.csv';
    end

    raw = csvread(filename, 1, 0);  % skip header row
    if isempty(raw)
        error('load_secondary_sync_dump: no data rows found in %s', filename);
    end

    scan_id_col = raw(:, 1);
    unique_ids = unique(scan_id_col);
    n = numel(unique_ids);

    stalls = struct('scan_id', {}, 'prn', {}, 'system_char', {}, 'channel', {}, ...
                     'secondary_code_length', {}, 'corr_value', {});

    for i = 1:n
        rows = raw(scan_id_col == unique_ids(i), :);
        [~, order] = sort(rows(:, 6));  % call_index, defensive (writer already emits in call order)
        rows = rows(order, :);

        stalls(i).scan_id                = rows(1, 1);
        stalls(i).prn                    = rows(1, 2);
        stalls(i).system_char            = rows(1, 3);
        stalls(i).channel                = rows(1, 4);
        stalls(i).secondary_code_length  = rows(1, 5);
        stalls(i).corr_value             = rows(:, 7)';
    end

    if nargout == 0
        printf('%6s %5s %4s %9s %10s\n', 'scan', 'PRN', 'ch', 'lockLen', 'nCalls');
        for i = 1:n
            printf('%6d %5d %4d %9d %10d\n', stalls(i).scan_id, stalls(i).prn, stalls(i).channel, ...
                stalls(i).secondary_code_length, numel(stalls(i).corr_value));
        end
        printf('\n%d stalls loaded from %s\n', n, filename);
        printf('PRNs present: %s\n', mat2str(unique([stalls.prn])));
        clear stalls;  % avoid echoing the struct array to the console when called as a script
    end
end
