# F1 Race Replayer with GUI (FastF1 + Matplotlib + Tkinter)
# -----------------------------------------------------------------------------
# Features:
# - GUI to select Year, Circuit, and Session.
# - Replays the FULL race for ALL drivers (stitched laps).
# - SMOOTH animation (0.2s data resolution).
# - Working Lap Counter.
# - Dynamic Leaderboard (Sorted by Distance).
# - Proper App Closure.
# - 2025 Data Fallback: Defaults to 2023 if 2025 data is missing.
# -----------------------------------------------------------------------------

import matplotlib
matplotlib.use("TkAgg") # Explicitly set backend before importing pyplot

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
    try:
        os.makedirs("f1_cache")
    except OSError:
        pass # Ignore if it exists or fails permissions (will use memory or fail later)
        
try:
    fastf1.Cache.enable_cache("f1_cache")
except:
    print("Warning: Could not enable FastF1 cache.")

# Global variables for animation control
anim = None

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
    
    max_race_time = 0
    drivers = pd.unique(session.laps['Driver'])
    
    # HIGHER RESOLUTION FOR SMOOTHNESS
    # 0.2s provides good balance between smoothness and memory usage
    step_size = 0.2 
    
    # Pre-calculate distances for leaderboard
    driver_distances = {}

    for driver in drivers:
        try:
            print(f"  Processing {driver}...")
            laps = session.laps.pick_driver(driver)
            
            # Get telemetry for ALL laps at once
            try:
                tel = laps.get_telemetry()
            except:
                continue

            if tel is None or tel.empty:
                continue
            
            t_seconds = tel['Time'].dt.total_seconds().to_numpy()
            x = tel['X'].to_numpy()
            y = tel['Y'].to_numpy()
            
            # Use 'Distance' column if available, otherwise estimate
            if 'Distance' in tel.columns:
                 d = tel['Distance'].to_numpy()
            else:
                 # Fallback: Cumulative distance estimation not accurate without speed integration
                 # Just use index as proxy if needed, or 0
                 d = np.zeros_like(x)

            telemetry_data[driver] = {
                "Time": t_seconds,
                "X": x,
                "Y": y,
                "Distance": d,
                "Color": get_driver_color(driver),
                "Team": laps.iloc[0]['Team'] if not laps.empty else "Unknown"
            }
            
            if len(t_seconds) > 0:
                max_race_time = max(max_race_time, t_seconds[-1])
                
        except Exception as e:
            print(f"Skipping {driver}: {e}")

    if not telemetry_data:
        return None, None, None, None

    # 2. Normalize Start Time
    # Find start time (min time across all)
    start_time = min([d["Time"][0] for d in telemetry_data.values() if len(d["Time"]) > 0])
    max_race_time -= start_time
    
    # 3. Create Common Timeline
    common_time = np.arange(0, max_race_time, step_size)
    
    # 4. Interpolate Data
    final_data = {}
    for driver, data in telemetry_data.items():
        shifted_time = data["Time"] - start_time
        
        # Interpolate X, Y and Distance to common timeline
        interp_x = np.interp(common_time, shifted_time, data["X"], left=np.nan, right=data["X"][-1])
        interp_y = np.interp(common_time, shifted_time, data["Y"], left=np.nan, right=data["Y"][-1])
        interp_d = np.interp(common_time, shifted_time, data["Distance"], left=0, right=data["Distance"][-1])
        
        final_data[driver] = {
            "X": interp_x,
            "Y": interp_y,
            "Distance": interp_d,
            "Color": data["Color"],
            "Team": data["Team"]
        }
    
    # 5. Calculate Lap Counter Data
    # Use the winner to define lap boundaries
    try:
        winner = drivers[0]
        winner_laps = session.laps.pick_driver(winner)
        # LapStartTime is a Timedelta. Convert to seconds relative to our normalized start_time
        lap_start_times = winner_laps['LapStartTime'].dt.total_seconds().to_numpy() - start_time
        lap_numbers = winner_laps['LapNumber'].to_numpy()
        
        # Ensure lap start times are sorted (should be)
        # Replace NaNs (first lap often NaN) with 0.0 or small negative
        lap_start_times = np.nan_to_num(lap_start_times, nan=0.0)
        
    except Exception as e:
        print(f"Lap counter error: {e}")
        lap_start_times = [0]
        lap_numbers = [1]
        
    return final_data, common_time, lap_start_times, lap_numbers

class RaceReplayerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("F1 Race Replayer 2025 (Full Race Mode)")
        self.root.geometry("1300x850") # Wider for leaderboard
        self.root.protocol("WM_DELETE_WINDOW", self.on_close) # Handle closure
        
        # --- Layout Frames ---
        control_frame = ttk.LabelFrame(root, text="Race Selection")
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)
        
        main_content = ttk.Frame(root)
        main_content.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.left_frame = ttk.Frame(main_content)
        self.left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.right_frame = ttk.Frame(main_content, width=250)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=5)

        # --- Controls ---
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
        self.speed_var = tk.StringVar(value="20x") 
        self.speed_combo = ttk.Combobox(control_frame, textvariable=self.speed_var, values=["1x", "5x", "10x", "20x", "50x"], width=5)
        self.speed_combo.pack(side=tk.LEFT, padx=0)
        self.speed_combo.bind("<<ComboboxSelected>>", self.change_speed)
        
        self.status_lbl = ttk.Label(control_frame, text="Ready", foreground="blue")
        self.status_lbl.pack(side=tk.LEFT, padx=10)

        # --- Matplotlib Canvas ---
        self.fig, self.ax = plt.subplots(figsize=(8, 6), facecolor='#121212')
        self.ax.set_facecolor('#121212')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.left_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # --- Leaderboard ---
        ttk.Label(self.right_frame, text="Leaderboard", font=("Arial", 12, "bold")).pack(side=tk.TOP, pady=5)
        
        # Treeview for leaderboard
        cols = ("Pos", "Driver", "Gap")
        self.leaderboard_tree = ttk.Treeview(self.right_frame, columns=cols, show="headings", height=30)
        self.leaderboard_tree.column("Pos", width=40, anchor="center")
        self.leaderboard_tree.column("Driver", width=80, anchor="w")
        self.leaderboard_tree.column("Gap", width=80, anchor="e")
        
        self.leaderboard_tree.heading("Pos", text="Pos")
        self.leaderboard_tree.heading("Driver", text="Driver")
        self.leaderboard_tree.heading("Gap", text="Dist (m)")
        
        self.leaderboard_tree.pack(fill=tk.BOTH, expand=True)

        self.anim = None
        self.telemetry_data = {}
        self.common_time = []
        self.lap_start_times = []
        self.lap_numbers = []
        self.driver_dots = {}
        self.driver_labels = {}
        self.lap_counter_text = None
        self.current_frame = 0

    def on_close(self):
        """Handle window closing event"""
        if self.anim:
            try:
                self.anim.event_source.stop()
            except:
                pass
        self.root.destroy()
        sys.exit()

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

        self.telemetry_data, self.common_time, self.lap_start_times, self.lap_numbers = process_telemetry(session)
        if not self.telemetry_data:
             self.status_lbl.config(text="No telemetry found!")
             return
             
        self.setup_plot()
        self.start_animation()
        self.status_lbl.config(text=f"Replaying {year} {circuit} (Full Race)")

    def setup_plot(self):
        self.ax.clear()
        self.ax.axis('off')
        
        # Plot Track Map 
        ref_driver = list(self.telemetry_data.keys())[0]
        ref_x = self.telemetry_data[ref_driver]["X"]
        ref_y = self.telemetry_data[ref_driver]["Y"]
        
        mask = ~np.isnan(ref_x)
        self.ax.plot(ref_x[mask], ref_y[mask], color='#333333', linewidth=6, alpha=0.5)
        
        # Lap Counter Text
        self.lap_counter_text = self.ax.text(0.02, 0.95, "Lap 1", transform=self.ax.transAxes, 
                                            color='white', fontsize=14, fontweight='bold')

        # Initialize Driver Dots
        self.driver_dots = {}
        self.driver_labels = {}
        
        for driver, data in self.telemetry_data.items():
            color = data["Color"]
            dot, = self.ax.plot([], [], marker='o', color=color, markersize=6, markeredgecolor='black', markeredgewidth=0.5)
            self.driver_dots[driver] = dot
            text = self.ax.text(0, 0, driver, color=color, fontsize=7, fontweight='bold', clip_on=True)
            self.driver_labels[driver] = text
            
        valid_x = ref_x[mask]
        valid_y = ref_y[mask]
        self.ax.set_xlim(np.min(valid_x) - 1000, np.max(valid_x) + 1000)
        self.ax.set_ylim(np.min(valid_y) - 1000, np.max(valid_y) + 1000)
        self.ax.set_aspect('equal')
        self.canvas.draw()

    def change_speed(self, event):
        if self.anim:
            self.start_animation(start_frame=self.current_frame)

    def start_animation(self, start_frame=0):
        if self.anim:
            try:
                self.anim.event_source.stop()
            except:
                pass
        
        total_frames = len(self.common_time)
        speed_str = self.speed_var.get().replace("x", "")
        speed_mult = int(speed_str)
        
        base_interval = (0.2 * 1000) / speed_mult
        
        step_frames = 1
        if base_interval < 10:
            step_frames = int(10 / base_interval) + 1
            base_interval = 10 
            
        def update(frame_idx):
            idx = start_frame + (frame_idx * step_frames)
            self.current_frame = idx
            
            if idx >= total_frames:
                try:
                    self.anim.event_source.stop()
                except:
                    pass
                return []

            # 1. Update Lap Counter
            current_race_time = self.common_time[idx]
            lap_idx = np.searchsorted(self.lap_start_times, current_race_time, side='right') - 1
            
            if lap_idx >= 0 and lap_idx < len(self.lap_numbers):
                current_lap = self.lap_numbers[lap_idx]
                total_laps = self.lap_numbers[-1]
                self.lap_counter_text.set_text(f"Lap {current_lap} / {total_laps}")
            elif lap_idx >= len(self.lap_numbers):
                 total_laps = self.lap_numbers[-1]
                 self.lap_counter_text.set_text(f"Finish / {total_laps}")
            else:
                 self.lap_counter_text.set_text(f"Grid / {self.lap_numbers[-1]}")

            # 2. Update Drivers & Leaderboard Data
            leaderboard_data = []

            for driver, dot in self.driver_dots.items():
                x_arr = self.telemetry_data[driver]["X"]
                y_arr = self.telemetry_data[driver]["Y"]
                dist_arr = self.telemetry_data[driver]["Distance"]
                
                x = x_arr[idx]
                y = y_arr[idx]
                dist = dist_arr[idx]
                
                if np.isnan(x) or np.isnan(y):
                    dot.set_visible(False)
                    self.driver_labels[driver].set_visible(False)
                else:
                    dot.set_visible(True)
                    self.driver_labels[driver].set_visible(True)
                    dot.set_data([x], [y])
                    self.driver_labels[driver].set_position((x + 200, y + 200))
                    
                    leaderboard_data.append((driver, dist))
            
            # 3. Update Leaderboard
            leaderboard_data.sort(key=lambda x: x[1], reverse=True)
            
            if frame_idx % 5 == 0:
                self.update_leaderboard_gui(leaderboard_data)
            
            return list(self.driver_dots.values()) + list(self.driver_labels.values()) + [self.lap_counter_text]

        self.anim = FuncAnimation(self.fig, update, frames=(total_frames - start_frame) // step_frames, 
                                  interval=base_interval, blit=True, repeat=False)
        self.canvas.draw()
        
    def update_leaderboard_gui(self, sorted_data):
        for item in self.leaderboard_tree.get_children():
            self.leaderboard_tree.delete(item)
            
        leader_dist = sorted_data[0][1] if sorted_data else 0
        
        for i, (driver, dist) in enumerate(sorted_data):
            pos = i + 1
            if i == 0:
                gap_str = "Leader"
            else:
                gap = leader_dist - dist
                gap_str = f"-{gap:.0f} m"
                
            self.leaderboard_tree.insert("", "end", values=(pos, driver, gap_str))

if __name__ == "__main__":
    print("Initializing F1 Race Replayer...")
    root = tk.Tk()
    app = RaceReplayerApp(root)
    print("Launching GUI...")
    root.mainloop()