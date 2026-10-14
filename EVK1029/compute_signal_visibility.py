#!/usr/bin/env python3
"""Compute GPS L1/L5 and Galileo E1/E5a signal visibility from a merged
RINEX 4.02 broadcast navigation file, at the AGNSS reference location used
by the EVK1029 configs, over a UTC time window.

Per-signal health is honored, not just per-satellite health: GPS carries a
distinct L1/L2/L5 health triplet (LNAV vs CNAV messages), and Galileo
carries a distinct E1B/E5a/E5b DVS+HS bitfield (INAV vs FNAV messages) --
a satellite can be healthy on one signal and not on another. A signal
whose health bits are not OK at a given epoch is excluded from that
signal's visibility, matching how gnss-sdr itself would refuse to use it.

Field layout and health bit meanings are taken directly from the RINEX
4.02 spec (IGS/RTCM RINEX Committee, 1 September 2025), Tables A9 (GPS
LNAV), A10 (GPS CNAV), A13 (Galileo INAV/FNAV).

Satellite position is always computed from the legacy message (LNAV for
GPS, INAV for Galileo) -- the broadcast orbit is the same physical orbit
regardless of which message carries it, so there is no need to trust
CNAV's/FNAV's own (differently-parameterized) orbit terms just to get a
per-signal health flag out of them.
"""
import argparse
import csv
import math
from collections import defaultdict
from datetime import datetime, timedelta, timezone

import matplotlib
matplotlib.use("Agg")  # headless -- always save to a file, never plt.show()
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Settings -- edit these and rerun, or override any of them on the command
# line (see --help).
# ---------------------------------------------------------------------------

# Satellites below this elevation are never counted as visible, regardless
# of health. gnss-sdr's own runtime-visibility-aware search and
# satellite_visibility.cc use the same El > elevation_mask_deg convention.
ELEVATION_MASK_DEG = 15.0

NAV_FILE = "/media/labusr/06F06D6F0CFF0EEC/evk1029_archive/BRD400DLR_S_20262530000_01D_MN.rnx"
DATE = "2026-09-10"
START_UTC = "07:27:00"
END_UTC = "10:00:00"
STEP_S = 60

# PRNs to exclude entirely, same semantics as gnss-sdr's own
# GNSS-SDR.Galileo_banned_prns / GNSS-SDR.GPS_banned_prns config keys
# (gnss_flowgraph.cc): removed from consideration outright, not just
# flagged unhealthy -- so this reference never counts a banned PRN as
# visible even if its broadcast health says OK. Comma-separated string,
# copy-pasteable straight from the .conf file's value (e.g.
# "1,14,17,20,22,24,35,32"); leave empty for none.
GALILEO_BANNED_PRNS = ""
GPS_BANNED_PRNS = ""

# AGNSS reference location used by the EVK1029 configs.
OBSERVER_LAT_DEG = 46.012965
OBSERVER_LON_DEG = 8.907565
OBSERVER_HEIGHT_M = 400.0

LEAP_SECONDS = 18  # GPST - UTC, per the nav file header's LEAP SECONDS line

OUTPUT_CSV = "signal_visibility_20260910_0727_1000.csv"
OUTPUT_PLOT = "signal_visibility_20260910_0727_1000.png"

# ---------------------------------------------------------------------------
# RINEX 4.02 parsing
# ---------------------------------------------------------------------------

# Orbit-1..7 field counts per nav type (how many of the 4-values-per-line
# broadcast orbit lines to read after the epoch/clock line).
ORBIT_LINES = {
    "LNAV": 7,
    "CNAV": 8,
    "INAV": 7,
    "FNAV": 7,
}


def _parse_rinex_float(tok):
    # RINEX floats use D or E exponent markers; D is not valid Python float
    # syntax.
    return float(tok.replace("D", "E").replace("d", "e"))


def _split_fixed(line, width=19):
    # Broadcast orbit lines are 4X followed by up to 4 fields of 19 chars
    # each (E19.12) -- fixed-width, not whitespace-delimited, since fields
    # can butt up against each other with no space (e.g. two negative
    # numbers in a row).
    body = line[4:]
    vals = []
    for i in range(0, len(body), width):
        chunk = body[i:i + width].strip()
        if chunk:
            vals.append(_parse_rinex_float(chunk))
    return vals


def parse_rinex_nav(path, wanted_types):
    """Yields dicts, one per EPH record, for nav message types in
    wanted_types (e.g. {"LNAV", "CNAV"} for GPS, {"INAV", "FNAV"} for
    Galileo)."""
    with open(path, "r") as f:
        lines = f.readlines()

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if line.startswith("> EPH"):
            parts = line.split()
            # "> EPH G01 LNAV" -> sat="G01", navtype="LNAV"
            sat = parts[2]
            navtype = parts[3]
            i += 1
            if navtype not in wanted_types or i >= n:
                # Skip this record's data lines (epoch line + orbit lines)
                # by scanning forward to the next "> " header.
                while i < n and not lines[i].startswith("> "):
                    i += 1
                continue

            epoch_line = lines[i]
            i += 1
            # The date fields ("G01 2026 09 10 00 00 00") are reliably
            # whitespace-separated, but the 3 clock fields right after them
            # are fixed-width (E19.12) with no guaranteed separator --
            # two adjacent negative values butt up against each other with
            # no space, same issue as the orbit lines below.
            date_part, clock_part = epoch_line[:23], epoch_line[23:]
            ep = date_part.split()
            year, month, day, hour, minute = (int(x) for x in ep[1:6])
            second = float(ep[6])
            af0, af1, af2 = (
                _parse_rinex_float(clock_part[j:j + 19])
                for j in (0, 19, 38)
            )

            n_orbit_lines = ORBIT_LINES[navtype]
            orbit = []
            for _ in range(n_orbit_lines):
                orbit.append(_split_fixed(lines[i]))
                i += 1
            yield {
                "sat": sat,
                "navtype": navtype,
                "epoch": datetime(year, month, day, hour, minute, int(second),
                                   tzinfo=timezone.utc),
                "af0": af0, "af1": af1, "af2": af2,
                "orbit": orbit,
            }
        else:
            i += 1


# ---------------------------------------------------------------------------
# Time helpers
# ---------------------------------------------------------------------------

GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)


def utc_to_gpst_seconds(dt_utc, leap_seconds):
    """Seconds since the GPS epoch (continuous, no week wrap) for a UTC
    datetime."""
    return (dt_utc - GPS_EPOCH).total_seconds() + leap_seconds


# ---------------------------------------------------------------------------
# Broadcast orbit -> ECEF (standard GPS ICD-200 / Galileo OS-SIS-ICD
# Keplerian propagation with second-harmonic perturbation terms -- the same
# algorithm, just different mu depending on constellation).
# ---------------------------------------------------------------------------

MU_GPS = 3.986005e14
MU_GAL = 3.986004418e14
OMEGA_E_DOT = 7.2921151467e-5


def eph_to_ecef(orbit, t_sow, mu):
    # orbit[0] = [IODE, Crs, dn0, M0]
    # orbit[1] = [Cuc, e, Cus, sqrtA]
    # orbit[2] = [Toe, Cic, OMEGA0, Cis]  -- Toe is the authoritative
    #            broadcast Time of Ephemeris (sec of week), used directly
    #            rather than re-derived from the epoch/Toc line.
    # orbit[3] = [i0, Crc, omega, OMEGADOT]
    # orbit[4] = [IDOT, ...]
    _, crs, dn0, m0 = orbit[0]
    cuc, e, cus, sqrtA = orbit[1]
    toe_sow, cic, omega0, cis = orbit[2]
    i0, crc, omega, omegadot = orbit[3]
    idot = orbit[4][0]

    A = sqrtA * sqrtA
    n0 = math.sqrt(mu / A ** 3)
    tk = t_sow - toe_sow
    if tk > 302400:
        tk -= 604800
    elif tk < -302400:
        tk += 604800

    n = n0 + dn0
    M = m0 + n * tk
    E = M
    for _ in range(12):
        E = E - (E - e * math.sin(E) - M) / (1 - e * math.cos(E))

    nu = math.atan2(math.sqrt(1 - e * e) * math.sin(E), math.cos(E) - e)
    phi = nu + omega
    du = cus * math.sin(2 * phi) + cuc * math.cos(2 * phi)
    dr = crs * math.sin(2 * phi) + crc * math.cos(2 * phi)
    di = cis * math.sin(2 * phi) + cic * math.cos(2 * phi)

    u = phi + du
    r = A * (1 - e * math.cos(E)) + dr
    inc = i0 + di + idot * tk
    Omega = omega0 + (omegadot - OMEGA_E_DOT) * tk - OMEGA_E_DOT * toe_sow

    xp = r * math.cos(u)
    yp = r * math.sin(u)

    X = xp * math.cos(Omega) - yp * math.cos(inc) * math.sin(Omega)
    Y = xp * math.sin(Omega) + yp * math.cos(inc) * math.cos(Omega)
    Z = yp * math.sin(inc)
    return X, Y, Z


# ---------------------------------------------------------------------------
# Geodetic <-> ECEF, az/el
# ---------------------------------------------------------------------------

WGS84_A = 6378137.0
WGS84_F = 1 / 298.257223563
WGS84_E2 = WGS84_F * (2 - WGS84_F)


def llh_to_ecef(lat_deg, lon_deg, h_m):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sin_lat = math.sin(lat)
    N = WGS84_A / math.sqrt(1 - WGS84_E2 * sin_lat * sin_lat)
    x = (N + h_m) * math.cos(lat) * math.cos(lon)
    y = (N + h_m) * math.cos(lat) * math.sin(lon)
    z = (N * (1 - WGS84_E2) + h_m) * sin_lat
    return x, y, z


def az_el(observer_llh, observer_ecef, sat_ecef):
    lat_deg, lon_deg, _ = observer_llh
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)

    dx = sat_ecef[0] - observer_ecef[0]
    dy = sat_ecef[1] - observer_ecef[1]
    dz = sat_ecef[2] - observer_ecef[2]
    rng = math.sqrt(dx * dx + dy * dy + dz * dz)

    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    sin_lon, cos_lon = math.sin(lon), math.cos(lon)

    east = -sin_lon * dx + cos_lon * dy
    north = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz
    up = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz

    az = math.degrees(math.atan2(east, north))
    if az < 0:
        az += 360.0
    el = math.degrees(math.asin(up / rng))
    return az, el


# ---------------------------------------------------------------------------
# Health decoding (RINEX 4.02 Table A9/A10/A13)
# ---------------------------------------------------------------------------

def gps_lnav_l1_healthy(orbit):
    svh = int(round(orbit[5][1]))
    return svh == 0


def gps_cnav_l5_healthy(orbit):
    # Matches gnss-sdr's own handling (rtklib_pvt_gs.cc: "if
    # (gps_cnav_ephemeris->signal_health != 0)") -- the whole 3-bit
    # L1/L2/L5 health value gates the signal, not just the L5 bit.
    svh = int(round(orbit[5][1]))
    return svh == 0


def gal_inav_e1_healthy(orbit):
    # Bit 0: E1B DVS, bits 1-2: E1B HS.
    svh = int(round(orbit[5][1]))
    return (svh & 0b0000111) == 0


def gal_fnav_e5a_healthy(orbit):
    # Bit 3: E5a DVS, bits 4-5: E5a HS.
    svh = int(round(orbit[5][1]))
    return ((svh >> 3) & 0b111) == 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_banned_prns(s):
    """Parse a comma-separated PRN list (same syntax as gnss-sdr's
    Galileo_banned_prns/GPS_banned_prns config keys, e.g.
    "1,14,17,20,22,24,35,32") into a set of ints. Tolerant enough to accept
    a value copy-pasted straight out of the .conf file, including a
    trailing ";..." note some of these lines carry (e.g.
    "...,35,32;18") -- everything from the first ";" onward is dropped
    before splitting. Non-numeric leftover tokens are skipped rather than
    raising, so a stray/trailing comma or comment fragment doesn't crash
    the script. Empty string -> empty set."""
    s = s.split(";", 1)[0]
    banned = set()
    for tok in s.split(","):
        tok = tok.strip()
        if tok.isdigit():
            banned.add(int(tok))
    return banned


def nearest_record(records, target_epoch, max_age_s=7200):
    best = None
    best_dt = None
    for rec in records:
        dt = abs((rec["epoch"] - target_epoch).total_seconds())
        if dt <= max_age_s and (best_dt is None or dt < best_dt):
            best, best_dt = rec, dt
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nav_file", nargs="?", default=NAV_FILE)
    parser.add_argument("--date", default=DATE)
    parser.add_argument("--start-utc", default=START_UTC)
    parser.add_argument("--end-utc", default=END_UTC)
    parser.add_argument("--step-s", type=int, default=STEP_S)
    parser.add_argument("--lat", type=float, default=OBSERVER_LAT_DEG)
    parser.add_argument("--lon", type=float, default=OBSERVER_LON_DEG)
    parser.add_argument("--height", type=float, default=OBSERVER_HEIGHT_M)
    parser.add_argument("--leap-seconds", type=int, default=LEAP_SECONDS)
    parser.add_argument("--elevation-mask", type=float, default=ELEVATION_MASK_DEG)
    parser.add_argument("--galileo-banned-prns", default=GALILEO_BANNED_PRNS,
                         help="Comma-separated Galileo PRNs to exclude entirely, e.g. "
                              "'1,14,17,20,22,24,35,32' -- same syntax as gnss-sdr's own "
                              "GNSS-SDR.Galileo_banned_prns")
    parser.add_argument("--gps-banned-prns", default=GPS_BANNED_PRNS,
                         help="Comma-separated GPS PRNs to exclude entirely -- same syntax "
                              "as gnss-sdr's own GNSS-SDR.GPS_banned_prns")
    parser.add_argument("-o", "--output", default=OUTPUT_CSV)
    parser.add_argument("--plot-output", default=OUTPUT_PLOT)
    args = parser.parse_args()
    elevation_mask = args.elevation_mask
    galileo_banned = parse_banned_prns(args.galileo_banned_prns)
    gps_banned = parse_banned_prns(args.gps_banned_prns)

    day = datetime.strptime(args.date, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    t_start = day + (datetime.strptime(args.start_utc, "%H:%M:%S") -
                      datetime(1900, 1, 1))
    t_end = day + (datetime.strptime(args.end_utc, "%H:%M:%S") -
                    datetime(1900, 1, 1))

    print(f"Parsing {args.nav_file} ...")
    gps_recs = defaultdict(lambda: defaultdict(list))    # navtype -> prn -> [rec]
    gal_recs = defaultdict(lambda: defaultdict(list))
    for rec in parse_rinex_nav(args.nav_file, {"LNAV", "CNAV"}):
        # QZSS also uses the "LNAV"/"CNAV" navtype labels -- filter to the
        # GPS ('G') satellite system specifically.
        if rec["sat"][0] == "G":
            gps_recs[rec["navtype"]][rec["sat"]].append(rec)
    for rec in parse_rinex_nav(args.nav_file, {"INAV", "FNAV"}):
        gal_recs[rec["navtype"]][rec["sat"]].append(rec)

    gps_prns = sorted(gps_recs["LNAV"].keys())
    gal_prns = sorted(gal_recs["INAV"].keys())

    # Same semantics as gnss-sdr's own Galileo_banned_prns/GPS_banned_prns:
    # removed from consideration outright, so a banned PRN never appears in
    # the output/plot at all (not merely marked unhealthy) -- matching what
    # gnss-sdr itself would actually search/use.
    if gps_banned:
        excluded = [p for p in gps_prns if int(p[1:]) in gps_banned]
        gps_prns = [p for p in gps_prns if int(p[1:]) not in gps_banned]
        print(f"GPS banned PRNs excluded: {excluded}")
    if galileo_banned:
        excluded = [p for p in gal_prns if int(p[1:]) in galileo_banned]
        gal_prns = [p for p in gal_prns if int(p[1:]) not in galileo_banned]
        print(f"Galileo banned PRNs excluded: {excluded}")

    print(f"GPS PRNs (LNAV): {len(gps_prns)}, Galileo PRNs (INAV): {len(gal_prns)}")

    observer_llh = (args.lat, args.lon, args.height)
    observer_ecef = llh_to_ecef(*observer_llh)

    rows = []
    t = t_start
    while t <= t_end:
        t_gpst = utc_to_gpst_seconds(t, args.leap_seconds)
        t_sow = t_gpst % 604800

        # GPS L1 (LNAV) and L5 (position from LNAV, health from CNAV)
        for prn in gps_prns:
            lnav = nearest_record(gps_recs["LNAV"][prn], t)
            if lnav is None:
                continue
            X, Y, Z = eph_to_ecef(lnav["orbit"], t_sow, MU_GPS)
            az, el = az_el(observer_llh, observer_ecef, (X, Y, Z))

            l1_ok = gps_lnav_l1_healthy(lnav["orbit"])
            rows.append((t, "GPS", prn, "L1", az, el, l1_ok, el > elevation_mask and l1_ok))

            cnav = nearest_record(gps_recs["CNAV"][prn], t)
            l5_ok = gps_cnav_l5_healthy(cnav["orbit"]) if cnav else False
            rows.append((t, "GPS", prn, "L5", az, el, l5_ok, el > elevation_mask and l5_ok))

        # Galileo E1 (INAV) and E5a (position from INAV, health from FNAV)
        for prn in gal_prns:
            inav = nearest_record(gal_recs["INAV"][prn], t)
            if inav is None:
                continue
            X, Y, Z = eph_to_ecef(inav["orbit"], t_sow, MU_GAL)
            az, el = az_el(observer_llh, observer_ecef, (X, Y, Z))

            e1_ok = gal_inav_e1_healthy(inav["orbit"])
            rows.append((t, "Galileo", prn, "E1", az, el, e1_ok, el > elevation_mask and e1_ok))

            fnav = nearest_record(gal_recs["FNAV"][prn], t)
            e5a_ok = gal_fnav_e5a_healthy(fnav["orbit"]) if fnav else False
            rows.append((t, "Galileo", prn, "E5a", az, el, e5a_ok, el > elevation_mask and e5a_ok))

        t += timedelta(seconds=args.step_s)

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["utc_time", "system", "prn", "signal", "azimuth_deg",
                    "elevation_deg", "health_ok", "visible"])
        for t, sys_, prn, sig, az, el, ok, vis in rows:
            w.writerow([t.strftime("%Y-%m-%d %H:%M:%S"), sys_, prn, sig,
                        f"{az:.2f}", f"{el:.2f}", int(ok), int(vis)])

    n_visible = sum(1 for r in rows if r[7])
    print(f"Wrote {len(rows)} rows ({n_visible} visible) to {args.output}")

    plot_visible_satellite_counts(rows, args.plot_output, elevation_mask)


def plot_visible_satellite_counts(rows, output_path, elevation_mask):
    """Number of visible satellites per epoch, one line per (system, signal),
    honoring the same elevation mask and per-signal health already baked
    into each row's 'visible' flag. Also plots a combined GPS+Galileo line:
    the number of distinct satellites (not signals) usable for a joint
    fix, one count per satellite via each constellation's primary signal
    (GPS L1, Galileo E1) rather than summing every per-signal line -- GPS
    L1+L5 and Galileo E1+E5a are the same physical satellite, so counting
    both would double-count it instead of reflecting geometry/DOP-relevant
    satellite availability."""
    series_order = [("GPS", "L1"), ("GPS", "L5"), ("Galileo", "E1"), ("Galileo", "E5a")]
    counts = {key: defaultdict(int) for key in series_order}
    for t, sys_, _prn, sig, _az, _el, _ok, vis in rows:
        if vis:
            counts[(sys_, sig)][t] += 1

    times = sorted({t for t, *_ in rows})

    # Combined GPS+Galileo: one count per satellite, from each
    # constellation's primary signal (GPS and Galileo PRNs never collide,
    # so GPS L1 count + Galileo E1 count is already the distinct-satellite
    # total, no set/union bookkeeping needed).
    combined_counts = {
        t: counts[("GPS", "L1")].get(t, 0) + counts[("Galileo", "E1")].get(t, 0)
        for t in times
    }

    fig, (ax_signals, ax_combined) = plt.subplots(
        2, 1, figsize=(12, 9), sharex=True,
        gridspec_kw={"height_ratios": [2, 1]})

    for key in series_order:
        sys_, sig = key
        y = [counts[key].get(t, 0) for t in times]
        ax_signals.plot(times, y, marker=".", label=f"{sys_} {sig}")

    ax_signals.set_ylabel("Number of visible satellites")
    ax_signals.set_title(f"Signal visibility (elevation mask = {elevation_mask:g} deg)")
    ax_signals.legend()
    ax_signals.grid(True, alpha=0.3)

    combined_y = [combined_counts[t] for t in times]
    ax_combined.plot(times, combined_y, marker="o", linewidth=2, color="black",
                      label="Combined GPS+Galileo")
    ax_combined.set_xlabel("UTC time")
    ax_combined.set_ylabel("Number of visible satellites")
    ax_combined.set_title("Combined GPS+Galileo (distinct satellites, primary signal)")
    ax_combined.legend()
    ax_combined.grid(True, alpha=0.3)

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    print(f"Wrote plot to {output_path}")


if __name__ == "__main__":
    main()
