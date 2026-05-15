# planetar-sat

Satellite imagery → AI boat detection → tracking microservice.

Fetches Sentinel-1 SAR scenes from Copernicus Data Space, runs a CFAR-based ship
detector over them, tracks detections across passes, and publishes results to
`planetar-broker` as native `zmesg` envelopes.

Built for the CH13 v2 proposal flagship demo (dark-vessel detection over the
Salish Sea) and for the `zdefence` v1 project.

## Layout

```
planetar-sat/
├── src/planetar_sat/
│   ├── fetch/        # Copernicus Data Space (CDSE) Sentinel-1 GRD fetcher
│   ├── detect/       # CFAR baseline SAR ship detector
│   ├── track/        # IoU multi-target tracker (persistent IDs)
│   ├── bus/          # zmesg encoder + TCP publisher to planetar-broker
│   ├── cli.py        # `planetar-sat fetch|detect|track|run`
│   └── config.py
├── tests/
└── docs/
```

## Quickstart

```bash
make install                            # pip install -e .[dev] into .venv
export CDSE_USERNAME=…                  # Copernicus Data Space credentials
export CDSE_PASSWORD=…

# live (last N days)
.venv/bin/planetar-sat run --aoi salish-sea --days 3

# historical date range
.venv/bin/planetar-sat run --aoi salish-sea --start 2024-08-15 --end 2024-08-22

# pure replay from cached scenes (no CDSE auth needed — great for demos)
.venv/bin/planetar-sat replay data/scenes/*.tif --acquired-at 2024-08-15T14:02:00Z
```

The `run` and `replay` subcommands chain fetch → detect → track → publish,
emitting THREE envelope types:

- `sar.chip` — per-detection machine envelope
- `track.update` — per-tracked-vessel machine envelope
- `chat.pac.sar-detections` — human-readable chat line, picked up by the
  `#sar-detections` channel in the running planetar-ui (which subscribes
  via `chatTopic()` → `chat.<server>.<channel>`)

Historical-replay mode stamps `created_at_ns` to the scene's acquisition
time so the UI renders messages at the right wall-clock — not as "now."

## Bus topics

| Topic         | Schema                                                                | Emitted by         |
|---------------|-----------------------------------------------------------------------|--------------------|
| `sar.chip`    | JSON: `{scene_id, lat, lon, snr, bbox_px, acquired_at_ns}`            | `detect/cfar.py`   |
| `track.update`| JSON: `{track_id, lat, lon, speed_kn, course_deg, last_seen_ns, ...}` | `track/tracker.py` |

Wire framing matches the broker: 4-byte little-endian length prefix, then a
`zmesg` envelope (see `src/planetar_sat/bus/zmesg.py` — pure-Python mirror of
`~/github/sness23/zmesg/zmesg.h`).

## Datasets in scope (CH13 1a)

- **Sentinel-1 GRD** (live ingress, free via CDSE) — **implemented**
- **xView3** — labeled-evaluation backbone (dark-vessel ground truth, Sentinel-1 + co-located AIS)
- **MarineCadastre AIS / Danish DMA / Global Fishing Watch** — historical AIS for `ais.gap` correlation
- **Alaska Satellite Facility** — US Sentinel-1 mirror + pre-Sentinel SAR archive (ERS-1/2, ALOS PALSAR)

See [`docs/DATA-SOURCES.md`](docs/DATA-SOURCES.md) for concrete URLs, auth
patterns, archive depth, and how each one lands on the bus.

See `~/github/planetarx/planetar/05-DATASETS.md` for the full modality list
across all five sensor lanes.

## Status

v0.1 — scaffold. CFAR baseline is honest (no fake numbers). Chip classifier
head is not in this commit; the proposal's `sar.chip` plan is CFAR-first +
classifier-port-from-xView3-reference. That port lives in a follow-on commit.
