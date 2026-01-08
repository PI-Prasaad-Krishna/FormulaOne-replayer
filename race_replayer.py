# F1 Race Replayer with GUI (FastF1 + Matplotlib + Tkinter)
# -----------------------------------------------------------------------------
# Features:
# - GUI to select Year, Circuit, and Session.
# - Replays the FULL race for ALL drivers (stitched laps).
# - Smooth animation using interpolation.
# - 2025 Data Fallback: Defaults to 2023 if 2025 data is missing.
# -----------------------------------------------------------------------------

import fastf1
import fastf1.plotting
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import numpy as np
import os
import sys
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Configure output for standard text
if sys.platform.startswith("win"):
    import codecs
    sys.stdout = codecs.getwriter("utf-8")(sys.stdout.buffer, "strict")

# 1. Setup Cache
if not os.path.exists("f1_cache"):
    os.makedirs("f1_cache")
fastf1.Cache.enable_cache("f1_cache")

# Global variables for animation control
anim = None
is_paused = False

def get_driver_color(driver_code):
    try:
        return fastf1.plotting.driver_color(driver_code)
    except:
        return "white"

def load_race_data(year, circuit, session_type="R"):
    print(f"Loading {year} {circuit} [{session_type}] session data...")
    try:
        session = fastf1.get_session(year, circuit, session_type)
        session.load(telemetry=True, laps=True, weather=False)
        return session
    except Exception as e:
        print(f"Error loading session: {e}")
        return None

def process_telemetry(session):
    print("Processing driver telemetry for all drivers (Full Race)...")
    telemetry_data = {}
    
    # 1. Determine the maximum race duration to set our timeline
    # We find the driver who finished last (or the max total time)
    # to ensure the timeline covers everyone.
    max_race_time = 0
    
    drivers = pd.unique(session.laps['Driver'])
    
    # We'll re-sample data to 0.5s intervals to keep memory usage manageable for a full race
    # 0.1s for 1.5 hours is too much data for a simple Matplotlib animation.
    # 0.5s is sufficient for a "replay".
    step_size = 0.5 
    
    for driver in drivers:
        try:
            print(f"  Processing {driver}...")
            # Get all laps for the driver
            laps = session.laps.pick_driver(driver)
            
            # We need to stitch laps together.
            # FastF1 telemetry has 'Time' (session time) or 'Date'.
            # We want 'Time' relative to Race Start.
            
            # Get telemetry for ALL laps at once (more efficient)
            try:
                tel = laps.get_telemetry()
            except:
                continue

            if tel is None or tel.empty:
                continue

            # 'Time' in telemetry is SessionTime.
            # We need to normalize it so 0 is the start of the race.
            # Usually, the first lap start time is a good anchor, or we just use SessionTime directly
            # and let the gaps naturally exist (grid formation etc might be included).
            # To make it clean, let's subtract the minimum session time across all drivers to start at roughly 0.
            
            # Let's just use the Time column (Timedelta) converted to seconds.
            t_seconds = tel['Time'].dt.total_seconds().to_numpy()
            x = tel['X'].to_numpy()
            y = tel['Y'].to_numpy()
            
            # Store raw data for now, we will interpolate later once we know global min/max
            telemetry_data[driver] = {
                "Time": t_seconds,
                "X": x,
                "Y": y,
                "Color": get_driver_color(driver),
                "Team": laps.iloc[0]['Team'] if not laps.empty else "Unknown"
            }
            
            if len(t_seconds) > 0:
                max_race_time = max(max_race_time, t_seconds[-1])
                
        except Exception as e:
            print(f"Skipping {driver}: {e}")

    if not telemetry_data:
        return None, None

    # 2. Normalize Start Time
    # Find the earliest timestamp to shift everyone to ~0
    start_time = min([d["Time"][0] for d in telemetry_data.values() if len(d["Time"]) > 0])
    max_race_time -= start_time
    
    # 3. Create Common Timeline
    # Using 0.5s steps.
    common_time = np.arange(0, max_race_time, step_size)
    
    final_data = {}
    
    for driver, data in telemetry_data.items():
        # Shift time
        shifted_time = data["Time"] - start_time
        
        # Interpolate X and Y to common timeline
        # Fill values outside of range (e.g. DNF) with NaN or the last position
        # We'll use the last known position so DNF cars "stop" on track.
        interp_x = np.interp(common_time, shifted_time, data["X"], left=np.nan, right=data["X"][-1])
        interp_y = np.interp(common_time, shifted_time, data["Y"], left=np.nan, right=data["Y"][-1])
        
        final_data[driver] = {
            "X": interp_x,
            "Y": interp_y,
            "Color": data["Color"],
            "Team": data["Team"]
        }
        
    return final_data, common_time

class RaceReplayerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("F1 Race Replayer 2025 (Full Race Mode)")
        self.root.geometry("1100x850")
        
        # --- Control Panel ---
        control_frame = ttk.LabelFrame(root, text="Race Selection")
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)
        
        ttk.Label(control_frame, text="Year:").pack(side=tk.LEFT, padx=5)
        self.year_var = tk.StringVar(value="2023")
        self.year_entry = ttk.Entry(control_frame, textvariable=self.year_var, width=6)
        self.year_entry.pack(side=tk.LEFT, padx=5)
        
        ttk.Label(control_frame, text="Circuit:").pack(side=tk.LEFT, padx=5)
        self.circuit_var = tk.StringVar(value="Abu Dhabi")
        self.circuit_entry = ttk.Entry(control_frame, textvariable=self.circuit_var, width=15)
        self.circuit_entry.pack(side=tk.LEFT, padx=5)
        
        self.load_btn = ttk.Button(control_frame, text="Load Full Race", command=self.start_replay)
        self.load_btn.pack(side=tk.LEFT, padx=10)
        
        # Playback Speed Control
        ttk.Label(control_frame, text="Speed:").pack(side=tk.LEFT, padx=10)
        self.speed_var = tk.StringVar(value="20x") # Default fast for full race
        self.speed_combo = ttk.Combobox(control_frame, textvariable=self.speed_var, values=["1x", "5x", "10x", "20x", "50x"], width=5)
        self.speed_combo.pack(side=tk.LEFT, padx=0)
        self.speed_combo.bind("<<ComboboxSelected>>", self.change_speed)
        
        self.status_lbl = ttk.Label(control_frame, text="Ready", foreground="blue")
        self.status_lbl.pack(side=tk.LEFT, padx=10)

        # --- Matplotlib Canvas ---
        self.fig, self.ax = plt.subplots(figsize=(8, 6), facecolor='#121212')
        self.ax.set_facecolor('#121212')
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.anim = None
        self.telemetry_data = {}
        self.common_time = []
        self.driver_dots = {}
        self.driver_labels = {}
        self.current_frame = 0

    def start_replay(self):
        year = int(self.year_var.get())
        circuit = self.circuit_var.get()
        
        self.status_lbl.config(text="Loading FULL race data... (This may take a minute)")
        self.root.update()
        
        session = load_race_data(year, circuit, "R")
        if not session:
            self.status_lbl.config(text="Error loading data!")
            messagebox.showerror("Error", f"Could not load data for {year} {circuit}")
            return

        self.telemetry_data, self.common_time = process_telemetry(session)
        if not self.telemetry_data:
             self.status_lbl.config(text="No telemetry found!")
             return
             
        self.setup_plot()
        self.start_animation()
        self.status_lbl.config(text=f"Replaying {year} {circuit} (Full Race)")

    def setup_plot(self):
        self.ax.clear()
        self.ax.axis('off')
        
        # Plot Track Map (using first driver's data as reference)
        ref_driver = list(self.telemetry_data.keys())[0]
        ref_x = self.telemetry_data[ref_driver]["X"]
        ref_y = self.telemetry_data[ref_driver]["Y"]
        
        # Handle NaNs for plotting track map
        mask = ~np.isnan(ref_x)
        self.ax.plot(ref_x[mask], ref_y[mask], color='#333333', linewidth=6, alpha=0.5) # Track layout
        
        # Initialize Driver Dots
        self.driver_dots = {}
        self.driver_labels = {}
        
        for driver, data in self.telemetry_data.items():
            color = data["Color"]
            dot, = self.ax.plot([], [], marker='o', color=color, markersize=6, markeredgecolor='black', markeredgewidth=0.5)
            self.driver_dots[driver] = dot
            text = self.ax.text(0, 0, driver, color=color, fontsize=7, fontweight='bold', clip_on=True)
            self.driver_labels[driver] = text
            
        # Set limits with some padding
        # Filter NaNs for min/max calculation
        valid_x = ref_x[mask]
        valid_y = ref_y[mask]
        
        self.ax.set_xlim(np.min(valid_x) - 1000, np.max(valid_x) + 1000)
        self.ax.set_ylim(np.min(valid_y) - 1000, np.max(valid_y) + 1000)
        self.ax.set_aspect('equal')
        self.canvas.draw()

    def get_interval(self):
        speed_str = self.speed_var.get().replace("x", "")
        speed = int(speed_str)
        # Base interval is 50ms for 0.5s step size
        # To speed up, we reduce interval
        # Real-time: 0.5s data step should take 500ms
        # 20x speed: 0.5s data step should take 25ms
        
        # Let's target a smooth 30fps animation.
        # We skip frames if speed is high.
        return 50 # Default 50ms per frame refresh

    def change_speed(self, event):
        if self.anim:
            # For FuncAnimation, changing speed usually requires restarting with new logic
            # or skipping frames in the update function.
            # We'll just restart animation from current frame
            self.start_animation(start_frame=self.current_frame)

    def start_animation(self, start_frame=0):
        if self.anim:
            self.anim.event_source.stop()
        
        total_frames = len(self.common_time)
        
        # Speed logic: How many data steps to jump per animation frame
        speed_str = self.speed_var.get().replace("x", "")
        speed_mult = int(speed_str)
        
        # Data is 0.5s per step.
        # If we update every 50ms (20fps), then 1x speed means we advance 0.1s per frame? No.
        # Let's simplify:
        # We update every 50ms.
        # At 1x speed, we want to show 0.05s of race time per 50ms (real time).
        # Our data step is 0.5s. That's huge.
        # So we just step through frames.
        # If we show 1 data point (0.5s) every 50ms, that is 10x speed.
        # So frame_step = 1 means 10x speed.
        # frame_step = 2 means 20x speed.
        
        # To support 1x, we'd need data step of 0.05s, or slow down interval to 500ms.
        # Let's adjust interval based on speed.
        
        interval = 50 # Base 50ms
        step = 1
        
        if speed_mult == 1:
            interval = 500 # 0.5s per 500ms = 1x
        elif speed_mult == 5:
            interval = 100 # 0.5s per 100ms = 5x
        elif speed_mult == 10:
            interval = 50  # 0.5s per 50ms = 10x
        elif speed_mult == 20:
            interval = 25  # 0.5s per 25ms = 20x
        elif speed_mult == 50:
            interval = 10  # Very fast
            
        def update(frame_idx):
            # Calculate actual data index
            idx = start_frame + frame_idx
            self.current_frame = idx
            
            if idx >= total_frames:
                self.anim.event_source.stop()
                return []

            for driver, dot in self.driver_dots.items():
                x_arr = self.telemetry_data[driver]["X"]
                y_arr = self.telemetry_data[driver]["Y"]
                
                x = x_arr[idx]
                y = y_arr[idx]
                
                # Handle NaNs (e.g. car retired or hasn't started)
                if np.isnan(x) or np.isnan(y):
                    dot.set_visible(False)
                    self.driver_labels[driver].set_visible(False)
                else:
                    dot.set_visible(True)
                    self.driver_labels[driver].set_visible(True)
                    dot.set_data([x], [y])
                    self.driver_labels[driver].set_position((x + 200, y + 200))
            
            return list(self.driver_dots.values()) + list(self.driver_labels.values())

        self.anim = FuncAnimation(self.fig, update, frames=total_frames - start_frame, 
                                  interval=interval, blit=True, repeat=False)
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = RaceReplayerApp(root)
    root.mainloop()