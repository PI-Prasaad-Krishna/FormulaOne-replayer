**F1 Visualizer Hub** 🏎️

An interactive Python tool suite for analyzing and visualizing Formula 1 races using telemetry from the FastF1 library.

**Overview**

This project consists of two main applications accessible via a central hub (`main.py`):
1.  **Race Replayer**: Replays "ghost cars" to visualize the race flow.
2.  **Analytics Dashboard**: A comprehensive suite for deep-dive analysis.

**Features**

### 📊 Analytics Dashboard (New!)
-   **Race Line Comparison**: Compare driver lines on a realistic track surface.
    -   **Apex Logic**: Visualizes the track ribbon with curvature-based shifting to show drivers hitting the apex.
    -   **Zoom Support**: Click to zoom in on corners for detailed analysis.
    -   **Smooth Visuals**: Uses B-spline interpolation for high-quality lines.
-   **Telemetry Battle**: Head-to-head speed and throttle/brake analysis.
-   **Race Strategy**: Visualizes lap time evolution and tyre history.
-   **Track Dominance**: Speed and gear shift maps.
-   **Driver DNA**: Radar charts comparing driver performance stats.
-   **Position Chart**: Lap-by-lap position changes.

### 🏁 Race Replayer
-   Full-race replay with smooth animation.
-   Adjustable playback speed (1x–50x).
-   Handles missing telemetry (DNFs).

**Requirements**

-   Python 3.8 or newer
-   Packages: `fastf1`, `matplotlib`, `pandas`, `customtkinter`, `scipy`

Install requirements:
```bash
pip install fastf1 matplotlib pandas customtkinter scipy
```

**Quick Start**

1.  Run the main hub:
```bash
python main.py
```
2.  Select **Race Replayer** or **Analytics Dashboard**.

**Acknowledgments**
Thanks to FastF1 and Matplotlib. Not affiliated with Formula 1.