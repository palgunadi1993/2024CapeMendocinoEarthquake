#!/usr/bin/env python3
"""Download a CMTSOLUTION file from the GlobalCMT NDK catalog.

Usage:
    download_cmtsolution.py --date 2024-12-05 --lat 40.35 --lon -125.0 --min-mag 6.5 -o data/cmtsolution

The script downloads the monthly NDK file from GlobalCMT, finds the closest
matching event by date, location, and magnitude, and writes it in CMTSOLUTION format.
"""

import argparse
import sys
import urllib.request
from datetime import datetime
from math import radians, cos, sin, asin, sqrt


GCMT_MONTHLY_URL = (
    "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_MONTHLY/{year}/{mon}{yr2}.ndk"
)
GCMT_QUICK_URL = (
    "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk"
)

MONTH_ABBR = [
    "jan", "feb", "mar", "apr", "may", "jun",
    "jul", "aug", "sep", "oct", "nov", "dec",
]


def haversine(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    return 2 * 6371 * asin(sqrt(a))


def download_ndk(year, month):
    """Download NDK catalog file. Try monthly first, fall back to quick CMT."""
    mon = MONTH_ABBR[month - 1]
    yr2 = f"{year % 100:02d}"
    url = GCMT_MONTHLY_URL.format(year=year, mon=mon, yr2=yr2)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            return resp.read().decode("utf-8")
    except Exception:
        pass
    # Fall back to quick CMT catalog
    print(f"Monthly catalog not found at {url}, trying quick CMT catalog...")
    with urllib.request.urlopen(GCMT_QUICK_URL, timeout=30) as resp:
        return resp.read().decode("utf-8")


def parse_ndk_events(ndk_text):
    """Parse NDK text into list of event dicts (5 lines per event)."""
    lines = ndk_text.strip().split("\n")
    events = []
    for i in range(0, len(lines) - 4, 5):
        block = lines[i : i + 5]
        if len(block) < 5:
            break
        try:
            events.append(parse_ndk_block(block))
        except Exception:
            continue
    return events


def parse_ndk_block(block):
    """Parse a 5-line NDK block into event dict."""
    h = block[0]
    src = h[:4].strip()
    parts = h[5:].split()
    yr, mo, dy = parts[0].split("/")
    time_parts = parts[1].split(":")
    hr = time_parts[0]
    mn = time_parts[1]
    sec = time_parts[2]
    lat = float(parts[2])
    lon = float(parts[3])
    dep = float(parts[4])
    mb = float(parts[5])
    ms = float(parts[6])
    region = " ".join(parts[7:])

    ev_name = block[1][:16].strip()
    hdur = float(block[1].split("TRIHD:")[1].strip())

    c = block[2].split()
    tshift = float(c[1])
    lat_c = float(c[3])
    lon_c = float(c[5])
    dep_c = float(c[7])

    m = block[3].split()
    exp = int(m[0])
    mrr = float(m[1])
    mtt = float(m[3])
    mpp = float(m[5])
    mrt = float(m[7])
    mrp = float(m[9])
    mtp = float(m[11])

    return {
        "src": src, "year": int(yr), "month": int(mo), "day": int(dy),
        "hour": int(hr), "minute": int(mn), "second": float(sec),
        "lat": lat, "lon": lon, "dep": dep, "mb": mb, "ms": ms,
        "region": region, "event_name": ev_name,
        "tshift": tshift, "hdur": hdur,
        "lat_c": lat_c, "lon_c": lon_c, "dep_c": dep_c,
        "exp": exp, "mrr": mrr, "mtt": mtt, "mpp": mpp,
        "mrt": mrt, "mrp": mrp, "mtp": mtp,
    }


def find_best_match(events, target_date, target_lat, target_lon, min_mag):
    """Find the event closest to the target parameters."""
    candidates = []
    for ev in events:
        mag = max(ev["mb"], ev["ms"])
        if mag < min_mag:
            continue
        ev_date = datetime(ev["year"], ev["month"], ev["day"])
        dt = abs((ev_date - target_date).total_seconds())
        dist = haversine(target_lat, target_lon, ev["lat"], ev["lon"])
        candidates.append((dt + dist * 100, ev))  # weight distance more

    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0])
    return candidates[0][1]


def write_cmtsolution(ev, outfile):
    """Write event in CMTSOLUTION format."""
    with open(outfile, "w") as f:
        # Source code and year must be joined (e.g. "PDEW2024") for
        # compatibility with neic-finitefault's read_gcmt_file parser.
        f.write(
            f" {ev['src']}{ev['year']:4d} {ev['month']:2d} {ev['day']:2d}"
            f" {ev['hour']:2d} {ev['minute']:2d}"
            f" {ev['second']:4.1f}"
            f" {ev['lat']:8.4f} {ev['lon']:9.4f} {ev['dep']:5.1f}"
            f" {ev['mb']:.1f} {ev['ms']:.1f}"
            f" {ev['region']}\n"
        )
        f.write(f"event name:     {ev['event_name']}\n")
        f.write(f"time shift:{ev['tshift']:11.4f}\n")
        f.write(f"half duration:{ev['hdur']:8.4f}\n")
        f.write(f"latitude:{ev['lat_c']:13.4f}\n")
        f.write(f"longitude:{ev['lon_c']:12.4f}\n")
        f.write(f"depth:{ev['dep_c']:16.4f}\n")
        e = ev["exp"]
        for name, val in [
            ("Mrr", ev["mrr"]), ("Mtt", ev["mtt"]), ("Mpp", ev["mpp"]),
            ("Mrt", ev["mrt"]), ("Mrp", ev["mrp"]), ("Mtp", ev["mtp"]),
        ]:
            f.write(f"{name}:{val:17.6f}e+{e:02d}\n")


def main():
    parser = argparse.ArgumentParser(description="Download CMTSOLUTION from GlobalCMT")
    parser.add_argument("--date", required=True, help="Event date (YYYY-MM-DD)")
    parser.add_argument("--lat", required=True, type=float, help="Approximate latitude")
    parser.add_argument("--lon", required=True, type=float, help="Approximate longitude")
    parser.add_argument("--min-mag", type=float, default=5.0, help="Minimum magnitude filter")
    parser.add_argument("-o", "--output", default="data/cmtsolution", help="Output file path")
    args = parser.parse_args()

    target_date = datetime.strptime(args.date, "%Y-%m-%d")

    print(f"Downloading GlobalCMT catalog for {target_date:%Y/%m}...")
    ndk_text = download_ndk(target_date.year, target_date.month)
    events = parse_ndk_events(ndk_text)
    print(f"Found {len(events)} events in catalog")

    ev = find_best_match(events, target_date, args.lat, args.lon, args.min_mag)
    if ev is None:
        print(f"No matching event found (date={args.date}, lat={args.lat}, lon={args.lon}, min_mag={args.min_mag})")
        sys.exit(1)

    print(f"Best match: {ev['event_name']} M{max(ev['mb'], ev['ms']):.1f} {ev['region']}")
    print(f"  Date: {ev['year']}-{ev['month']:02d}-{ev['day']:02d} {ev['hour']:02d}:{ev['minute']:02d}:{ev['second']:04.1f}")
    print(f"  Location: {ev['lat']:.2f}, {ev['lon']:.2f}, depth={ev['dep_c']:.1f} km")

    write_cmtsolution(ev, args.output)
    print(f"CMTSOLUTION written to {args.output}")


if __name__ == "__main__":
    main()
