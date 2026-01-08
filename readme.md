**F1 Telemetry Race Replayer** 🏎️

An interactive Python tool that visualizes Formula 1 races by replaying "ghost cars" using telemetry from the FastF1 library.

**Overview**

This project replays full race telemetry as a live, animated "ghost race" so you can observe gaps, cornering speed, and battles for position over time.

The app is optimized for full race sessions and works with cached session files stored in the `f1_cache/` folder (see the repository layout).

**Features**

- Interactive GUI for selecting year, event, and playback speed (tkinter).
- Full-race replay with smooth animation via interpolation.
- Adjustable playback speed (e.g., 1x–50x).
- Handles missing telemetry (DNFs) and normalizes start times.

**Requirements**

- Python 3.8 or newer
- Packages: `fastf1`, `matplotlib`, `pandas`
- `tkinter` (usually included with standard Python on Windows)

Install requirements with pip:

```bash
pip install fastf1 matplotlib pandas
```

**Quick Start**

1. Place or cache FastF1 session files under the `f1_cache/` folder (this repository already includes example cached sessions).
2. Run the replayer:

```bash
python race_replayer.py
```

3. In the GUI:
- Enter the **Year** (for example, `2023`).
- Enter the **Circuit** or event name (e.g., `Abu Dhabi`).
- Choose a **Speed** multiplier.
- Click **Load Full Race** to start the replay.

**Notes & Troubleshooting**

- If `tkinter` isn't available on your Python install, install the appropriate system package or use a Python installer that includes it.
- The app will use cached files in `f1_cache/` when available — this avoids repeated downloads via FastF1.

**How it works (brief)**

- Uses FastF1 to obtain telemetry and lap data.
- Extracts and normalizes telemetry per driver so all sessions align to a common timeline.
- Re-samples to a regular timeline and animates driver positions using Matplotlib's `FuncAnimation` embedded in a tkinter window.

**Future Improvements**

- Add a lap counter and live leaderboard.
- Smoother animation (soon ig? I'll try)
- Visualize pit stops and additional session types.
- Allow selecting subsets of drivers for focused replays.

**Acknowledgments**

This is a fan project — not affiliated with Formula 1, the FIA, or official FastF1 maintainers. Thanks to the FastF1 and Matplotlib projects for the tooling and data access.