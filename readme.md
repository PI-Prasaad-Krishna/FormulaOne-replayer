# F1 Visualizer Hub 🏎️

An interactive Python tool suite for analyzing and visualizing Formula 1 races using telemetry from the FastF1 library.

## Overview

This project provides a comprehensive set of tools for F1 data analysis, accessible via a central hub (`main.py`). It consists of two powerful applications:

1.  **Race Replayer**: A modern race visualization tool.
2.  **Analytics Dashboard**: A deep-dive telemetry and strategy analysis suite.

## Features

### 🏁 Race Replayer (Modern Edition)
-   **Full Race Animation**: Smooth 20-driver playback with interpolated positions.
-   **3D View**:
    -   **Soft Camera**: Physics-based navigation with inertia and smooth damping.
    -   **Ghost Curtain**: Subtle elevation visualization for track depth without visual clutter.
    -   **Optimized Rendering**: Smart blitting logic ensures silky smooth rotation even at 60FPS.
-   **Zero-Flicker Leaderboard**: High-performance UI with state caching for instant, jitter-free updates.
-   **Playback Control**: Pause, Resume, and adjustable speeds (1x–50x).
-   **Tyre Strategy**: Visual indicator of current tyre compounds for each driver.
-   **Session Info**: Displays real-time track status (SC, VSC, Red Flag) and weather conditions.
-   **Smart Movement Detection**: Automatically detects the race start to skip pre-race waiting times.

### 📊 Analytics Dashboard
-   **Telemetry Battle**: Head-to-head speed and throttle/brake trace comparison.
-   **Race Strategy**: Visualize lap time evolution and tyre history stints.
-   **Race Lines**: Compare driver racing lines on a realistic track surface with apex visualization.
-   **Tyre Degradation**: Analyze lap time drop-off vs. tyre age.
-   **Driver DNA**: Radar charts comparing drivers on Quali Pace, Race Pace, Consistency, and Tyre Management.
-   **Track Dominance**: Visual maps of speed and gear shifts around the circuit.
-   **Position Chart**: Lap-by-lap position changes for the entire grid.

## Installation

### 1. Requirements
Ensure you have Python 3.8+ installed.

### 2. Install Dependencies
Install the required Python packages using pip:

```bash
pip install fastf1 customtkinter matplotlib pandas scipy numpy
```

### 3. Setup
Clone the repository or download the source code:

```bash
git clone <repository-url>
cd "F1 Visulaizer"
```

> **Note**: On first run, a `f1_cache` directory will be created to store downloaded telemetry data.

## Quick Start

1.  Launch the main hub:
    ```bash
    python main.py
    ```
2.  Select your desired tool:
    -   Click **RACE REPLAYER** to watch a race unfold.
    -   Click **ANALYTICS DASHBOARD** to dive into data.

3.  **In the Apps**:
    -   Enter the **Year** (e.g., `2023`) and **Circuit** (e.g., `Abu Dhabi`).
    -   For comparisons, enter driver abbreviations (e.g., `VER` vs `HAM`).
    -   **3D Controls**: Drag to rotate (with inertia), Scroll to zoom.

## Troubleshooting

-   **High DPI Displays**: If the window looks too small or text is blurry, check your OS display scaling settings. CustomTkinter generally handles scaling well.
-   **Slow Loading**: First-time data loading for a race can take a minute as `fastf1` downloads and caches large telemetry files. Subsequent loads will be much faster.
-   **Missing Data**: Some older races or practice sessions might have incomplete telemetry.

## Acknowledgments
Powered by [FastF1](https://github.com/theOehrly/Fast-F1)
Not affiliated with Formula 1.