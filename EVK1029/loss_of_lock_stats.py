#!/usr/bin/env python3
"""Run gnss-sdr repeatedly against a config file and tally Loss-of-lock events by PRN.

Each trial launches gnss-sdr and lets it run until either the *receiver* (signal) time
printed on the console reaches --timeout seconds, or (if --min-signals-for-fix is given)
a position fix using at least that many satellites is seen -- whichever happens first,
however long that takes in wall clock time -- then kills it and parses whatever console
output was produced.
Receiver time comes from lines like:

    Current receiver time: 42 s
    Current receiver time: 1 min 3 s
    Current receiver time: 2 h 5 min 0 s

printed once per second of processed signal time (gnss_sdr_sample_counter), which is
NOT the same as wall-clock time since a file source with no throttle processes faster
or slower than real time depending on CPU load.

Within that captured output, lines like:

    Loss of lock in channel 3, satellite Galileo PRN E05 (Block I-A) !

are matched and counted per PRN, restricted to channels <= --max-channel. The last

    Position at ... using 7 observations is Lat = ...

line seen before the trial is killed gives the number of signals used for the fix
at timeout. The *first* such line seen in a trial also gives TTFF (time-to-first-fix),
reported in receiver time (seconds since the start of the run, from the same
"Current receiver time: ..." ticks described above) -- not wall-clock time, since
wall-clock TTFF would depend on CPU load / how fast the file source is read, not on
the receiver's own behavior.
"""

import argparse
import re
import statistics
import subprocess
import sys
from collections import Counter

ANSI_RE = re.compile(r"\x1B\[[0-?]*[ -/]*[@-~]")
LOSS_OF_LOCK_RE = re.compile(
    r"Loss of lock in channel (\d+),\s*satellite (\S+) PRN (\S+)"
)
POSITION_FIX_RE = re.compile(r"Position at .+? using (\d+) observations")
RX_TIME_RE = re.compile(
    r"Current receiver time:\s*"
    r"(?:(\d+)\s+days?\s+)?"
    r"(?:(\d+)\s+h\s+)?"
    r"(?:(\d+)\s+min\s+)?"
    r"(\d+)\s+s"
)


def parse_rx_time_seconds(line):
    """Return total receiver-time seconds from a 'Current receiver time: ...' line, or None."""
    match = RX_TIME_RE.search(line)
    if not match:
        return None
    days, hours, minutes, seconds = match.groups()
    total = int(seconds)
    total += int(minutes) * 60 if minutes else 0
    total += int(hours) * 3600 if hours else 0
    total += int(days) * 86400 if days else 0
    return total


def run_trial(exe, conf, workdir, rx_timeout_s, min_signals_for_fix=None):
    """Run one gnss-sdr trial until receiver time reaches rx_timeout_s, or (if
    min_signals_for_fix is given) a position fix using at least that many satellites is
    seen -- whichever comes first -- then kill it.

    Returns (output, stop_reason, ttff_rx_s, ttff_sv) where stop_reason is one of:
      "rx_time_limit"   - receiver time reached rx_timeout_s
      "fix_reached"     - a position fix with >= min_signals_for_fix satellites was seen
      "process_exited"  - gnss-sdr exited on its own before either condition was met
    ttff_rx_s is the receiver time (seconds, from the most recent "Current receiver
    time: ..." line seen so far) at which the first "Position at ..." line appeared, and
    ttff_sv is the number of observations used in that same first fix -- both None if no
    fix was seen during the trial.
    """
    cmd = [exe, "-c", conf]
    proc = subprocess.Popen(
        cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )

    stop_reason = "process_exited"
    last_rx_seconds = None
    ttff_rx_s = None
    ttff_sv = None

    lines = []
    try:
        for line in proc.stdout:
            lines.append(line)
            rx_seconds = parse_rx_time_seconds(line)
            if rx_seconds is not None:
                last_rx_seconds = rx_seconds
            if ttff_rx_s is None:
                first_fix_match = POSITION_FIX_RE.search(line)
                if first_fix_match is not None:
                    # First fix of the trial: receiver time at this point, from the most
                    # recent "Current receiver time: ..." tick (falls back to 0 if the fix
                    # somehow prints before the first tick, e.g. sub-second TTFF), plus how
                    # many observations that first fix used.
                    ttff_rx_s = last_rx_seconds if last_rx_seconds is not None else 0
                    ttff_sv = int(first_fix_match.group(1))
            if rx_seconds is not None and rx_seconds >= rx_timeout_s:
                stop_reason = "rx_time_limit"
                proc.kill()
                break
            if min_signals_for_fix is not None:
                fix_match = POSITION_FIX_RE.search(line)
                if fix_match is not None and int(fix_match.group(1)) >= min_signals_for_fix:
                    stop_reason = "fix_reached"
                    proc.kill()
                    break
    finally:
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()

    return ANSI_RE.sub("", "".join(lines)), stop_reason, ttff_rx_s, ttff_sv


def parse_loss_of_lock(output, max_channel):
    """Return list of (channel, prn) tuples for loss-of-lock events with channel <= max_channel."""
    events = []
    for match in LOSS_OF_LOCK_RE.finditer(output):
        channel = int(match.group(1))
        prn = match.group(3)
        if channel <= max_channel:
            events.append((channel, prn))
    return events


def parse_signals_at_timeout(output):
    """Return the number of observations used in the last position fix printed, or None if no fix was seen."""
    matches = POSITION_FIX_RE.findall(output)
    if not matches:
        return None
    return int(matches[-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conf", default="./EVK1029_Galileo_E1_E5a.conf", help="Config file passed to gnss-sdr -c")
    parser.add_argument("--exe", default="../install/gnss-sdr", help="Path to the gnss-sdr executable")
    parser.add_argument("--workdir", default=".", help="Working directory to run gnss-sdr from")
    parser.add_argument("--trials", type=int, default=50, help="Number of trials to run")
    parser.add_argument("--timeout", type=float, default=60.0, help="Receiver (signal) time in seconds, from 'Current receiver time:', to let each trial run before killing it -- always waited for, however long that takes")
    parser.add_argument("--max-channel", type=int, default=12, help="Only count loss-of-lock events on channels <= this value")
    parser.add_argument("--min-signals-for-fix", type=int, default=None, help="If set, stop a trial early (before --timeout) as soon as a position fix using at least this many satellites is seen")
    args = parser.parse_args()

    prn_counts = Counter()
    trials_with_events = 0
    total_events = 0
    signals_at_timeout = []  # one entry per trial that had a position fix
    ttff_values = []  # one entry per trial that had a position fix, in receiver-time seconds
    ttff_sv_values = []  # one entry per trial that had a position fix, # SV used at that first fix

    stop_reason_labels = {
        "rx_time_limit": "timeout",
        "fix_reached": f"SV=={args.min_signals_for_fix}",
        "process_exited": "process_exited",
    }

    for trial in range(1, args.trials + 1):
        output, stop_reason, ttff_rx_s, ttff_sv = run_trial(
            args.exe, args.conf, args.workdir, args.timeout, args.min_signals_for_fix
        )
        events = parse_loss_of_lock(output, args.max_channel)
        n_signals = parse_signals_at_timeout(output)

        if events:
            trials_with_events += 1
        total_events += len(events)
        for _channel, prn in events:
            prn_counts[prn] += 1

        if n_signals is not None:
            signals_at_timeout.append(n_signals)

        if ttff_rx_s is not None:
            ttff_values.append(ttff_rx_s)
        if ttff_sv is not None:
            ttff_sv_values.append(ttff_sv)

        prn_summary = ", ".join(f"{prn}" for _c, prn in events) if events else "none"
        signals_summary = str(n_signals) if n_signals is not None else "no fix"
        ttff_summary = f"{ttff_rx_s} s, {ttff_sv} SV" if ttff_rx_s is not None else "no fix"
        reason_summary = stop_reason_labels.get(stop_reason, stop_reason)
        print(
            f"Trial {trial}/{args.trials}: {len(events)} loss-of-lock event(s) [{prn_summary}], "
            f"TTFF: {ttff_summary}, stopped on: {reason_summary}, signals at stop: {signals_summary}",
            flush=True,
        )

    print()
    print(f"===== Statistics over {args.trials} trials (channels <= {args.max_channel}) =====")
    print(f"Trials with at least one loss-of-lock event: {trials_with_events}/{args.trials}")
    print(f"Total loss-of-lock events over the entire trials period: {total_events}")
    print()
    if prn_counts:
        print(f"{'PRN':<8}{'count':>8}{'% of trials':>14}")
        for prn, count in prn_counts.most_common():
            pct = 100.0 * count / args.trials
            print(f"{prn:<8}{count:>8}{pct:>13.1f}%")
    else:
        print("No loss-of-lock events observed for any PRN.")

    print()
    trials_with_fix = len(signals_at_timeout)
    print(f"Trials with a position fix at stop: {trials_with_fix}/{args.trials}")
    if signals_at_timeout:
        print(f"Signals at stop: mean={statistics.mean(signals_at_timeout):.2f}, "
              f"min={min(signals_at_timeout)}, max={max(signals_at_timeout)}, "
              f"median={statistics.median(signals_at_timeout):.1f}")
    else:
        print("No trial reached a position fix before stopping.")

    print()
    trials_with_ttff = len(ttff_values)
    print(f"Trials with TTFF (time-to-first-fix, receiver time): {trials_with_ttff}/{args.trials}")
    if ttff_values:
        print(f"TTFF: mean={statistics.mean(ttff_values):.2f} s, "
              f"min={min(ttff_values)} s, max={max(ttff_values)} s, "
              f"median={statistics.median(ttff_values):.1f} s")
        print(f"SV used at first fix: mean={statistics.mean(ttff_sv_values):.2f}, "
              f"min={min(ttff_sv_values)}, max={max(ttff_sv_values)}, "
              f"median={statistics.median(ttff_sv_values):.1f}")
    else:
        print("No trial reached a position fix before timeout.")


if __name__ == "__main__":
    sys.exit(main())
