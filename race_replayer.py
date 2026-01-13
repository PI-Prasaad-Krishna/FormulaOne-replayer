# F1 Race Replayer v2.0 (FastF1 + Matplotlib + CustomTkinter)
# -----------------------------------------------------------------------------
# Features:
# - MODERN GUI: Using CustomTkinter for a sleek Dark Mode look.
# - SMOOTH ANIMATION: 0.2s data resolution.
# - DYNAMIC LEADERBOARD: Updates in real-time without flickering.
# - FULL RACE SUPPORT: Stitches laps for all 20 drivers.
# -----------------------------------------------------------------------------

import matplotlib
matplotlib.use("TkAgg") # Required for integrating with Tkinter/CustomTkinter

import fastf1
import fastf1.plotting
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import numpy as np
import os
import sys
import customtkinter as ctk
import tkinter as tk
from tkinter import messagebox

# Configure output for standard text
if sys.platform.startswith("win"):
    import codecs
    sys.stdout = codecs.getwriter("utf-8")(sys.stdout.buffer, "strict")

# 1. Setup Cache
if not os.path.exists("f1_cache"):
    try:
        os.makedirs("f1_cache")
    except OSError:
        pass
        
try:
    fastf1.Cache.enable_cache("f1_cache")
except:
    print("Warning: Could not enable FastF1 cache.")

# 2. Configuration & Utilities
ctk.set_appearance_mode("Dark")  # Modes: "System" (standard), "Dark", "Light"
ctk.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"

def get_driver_color(driver_code):
    try:
        return fastf1.plotting.driver_color(driver_code)
    except:
        return "#ffffff"

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
    print("Processing driver telemetry...")
    telemetry_data = {}
    max_race_time = 0
    drivers = pd.unique(session.laps['Driver'])
    
    # High resolution for smoothness
    step_size = 0.2 
    
    for driver in drivers:
        try:
            laps = session.laps.pick_driver(driver)
            try:
                tel = laps.get_telemetry()
            except:
                continue

            if tel is None or tel.empty:
                continue
            
            t_seconds = tel['Time'].dt.total_seconds().to_numpy()
            x = tel['X'].to_numpy()
            y = tel['Y'].to_numpy()
            
            # Use Distance column if valid, else 0
            d = tel['Distance'].to_numpy() if 'Distance' in tel.columns else np.zeros_like(x)

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
                
        except Exception:
            pass

    if not telemetry_data:
        return None, None, None, None

    # Normalize Start Time
    start_time = min([d["Time"][0] for d in telemetry_data.values() if len(d["Time"]) > 0])
    max_race_time -= start_time
    
    # Create Common Timeline
    common_time = np.arange(0, max_race_time, step_size)
    
    # Interpolate Data
    final_data = {}
    for driver, data in telemetry_data.items():
        shifted_time = data["Time"] - start_time
        
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
    
    # Calculate Lap Data
    try:
        winner = drivers[0]
        winner_laps = session.laps.pick_driver(winner)
        lap_start_times = winner_laps['LapStartTime'].dt.total_seconds().to_numpy() - start_time
        lap_numbers = winner_laps['LapNumber'].to_numpy()
        lap_start_times = np.nan_to_num(lap_start_times, nan=0.0)
    except:
        lap_start_times = [0]
        lap_numbers = [1]
        
    return final_data, common_time, lap_start_times, lap_numbers

class RaceReplayerApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("F1 Race Replayer 2025 - Modern Edition")
        self.geometry("1400x900")
        self.protocol("WM_DELETE_WINDOW", self.on_close)

        # -- Layout Configuration --
        self.grid_columnconfigure(1, weight=1) # Main canvas expands
        self.grid_rowconfigure(1, weight=1)    # Main content expands

        # -- 1. Header / Controls --
        self.control_frame = ctk.CTkFrame(self, height=60, corner_radius=0)
        self.control_frame.grid(row=0, column=0, columnspan=3, sticky="ew", padx=0, pady=0)
        
        ctk.CTkLabel(self.control_frame, text="Year:", font=("Roboto", 14)).pack(side="left", padx=(20, 5))
        self.year_entry = ctk.CTkEntry(self.control_frame, width=60, placeholder_text="2023")
        self.year_entry.insert(0, "2023")
        self.year_entry.pack(side="left", padx=5)
        
        ctk.CTkLabel(self.control_frame, text="Circuit:", font=("Roboto", 14)).pack(side="left", padx=(20, 5))
        self.circuit_entry = ctk.CTkEntry(self.control_frame, width=150, placeholder_text="Abu Dhabi")
        self.circuit_entry.insert(0, "Abu Dhabi")
        self.circuit_entry.pack(side="left", padx=5)
        
        self.load_btn = ctk.CTkButton(self.control_frame, text="LOAD RACE", command=self.start_replay, fg_color="#E10600", font=("Roboto", 12, "bold"))
        self.load_btn.pack(side="left", padx=20)
        
        ctk.CTkLabel(self.control_frame, text="Speed:", font=("Roboto", 14)).pack(side="left", padx=(20, 5))
        self.speed_var = ctk.StringVar(value="20x")
        self.speed_combo = ctk.CTkComboBox(self.control_frame, values=["1x", "5x", "10x", "20x", "50x"], variable=self.speed_var, width=80, command=self.change_speed)
        self.speed_combo.pack(side="left", padx=5)
        
        self.status_lbl = ctk.CTkLabel(self.control_frame, text="Ready", text_color="gray")
        self.status_lbl.pack(side="left", padx=20)

        # -- 2. Leaderboard (Left Side) --
        self.leaderboard_frame = ctk.CTkFrame(self, width=250, corner_radius=0)
        self.leaderboard_frame.grid(row=1, column=0, sticky="nsew", padx=0, pady=0)
        
        ctk.CTkLabel(self.leaderboard_frame, text="LEADERBOARD", font=("Roboto", 16, "bold"), text_color="#E10600").pack(pady=(15,10))
        
        # Header Row
        header_row = ctk.CTkFrame(self.leaderboard_frame, fg_color="transparent")
        header_row.pack(fill="x", padx=10, pady=2)
        ctk.CTkLabel(header_row, text="POS", width=30, anchor="w", font=("Roboto", 10, "bold")).pack(side="left")
        ctk.CTkLabel(header_row, text="DRIVER", width=60, anchor="w", font=("Roboto", 10, "bold")).pack(side="left", padx=5)
        ctk.CTkLabel(header_row, text="GAP", width=60, anchor="e", font=("Roboto", 10, "bold")).pack(side="right")

        # Scrollable area for drivers
        self.scroll_lb = ctk.CTkScrollableFrame(self.leaderboard_frame, fg_color="transparent")
        self.scroll_lb.pack(fill="both", expand=True, padx=5, pady=5)
        
        # Pre-create 22 rows for performance (updating text is faster than creating widgets)
        self.lb_rows = []
        for i in range(22):
            row_frame = ctk.CTkFrame(self.scroll_lb, height=25, fg_color=("gray85", "gray20"))
            row_frame.pack(fill="x", pady=2)
            
            pos_lbl = ctk.CTkLabel(row_frame, text=f"{i+1}", width=30, font=("Roboto", 12, "bold"))
            pos_lbl.pack(side="left", padx=5)
            
            drv_lbl = ctk.CTkLabel(row_frame, text="-", font=("Roboto", 12), anchor="w")
            drv_lbl.pack(side="left", padx=5, fill="x", expand=True)
            
            gap_lbl = ctk.CTkLabel(row_frame, text="", font=("Mono", 11), width=60, anchor="e")
            gap_lbl.pack(side="right", padx=5)
            
            self.lb_rows.append((row_frame, pos_lbl, drv_lbl, gap_lbl))

        # -- 3. Map Canvas (Center/Right) --
        self.map_frame = ctk.CTkFrame(self, corner_radius=0, fg_color="#121212")
        self.map_frame.grid(row=1, column=1, columnspan=2, sticky="nsew")
        
        # Setup Matplotlib
        plt.style.use('dark_background')
        self.fig, self.ax = plt.subplots(figsize=(8, 6), facecolor='#121212')
        self.ax.set_facecolor('#121212')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.map_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # Animation vars
        self.anim = None
        self.telemetry_data = {}
        self.common_time = []
        self.driver_dots = {}
        self.driver_labels = {}
        self.lap_counter_text = None
        self.current_frame = 0

    def on_close(self):
        if self.anim:
            try: self.anim.event_source.stop()
            except: pass
        self.quit()
        self.destroy()
        sys.exit()

    def start_replay(self):
        try:
            year = int(self.year_entry.get())
            circuit = self.circuit_entry.get()
        except ValueError:
            messagebox.showerror("Input Error", "Year must be a number.")
            return
        
        self.status_lbl.configure(text="Loading Data... (Please Wait)", text_color="orange")
        self.update()
        
        session = load_race_data(year, circuit, "R")
        if not session:
            self.status_lbl.configure(text="Data Load Failed", text_color="red")
            messagebox.showerror("Error", f"Could not load data for {year} {circuit}")
            return

        self.telemetry_data, self.common_time, self.lap_start_times, self.lap_numbers = process_telemetry(session)
        if not self.telemetry_data:
             self.status_lbl.configure(text="No Telemetry Found", text_color="red")
             return
             
        self.setup_plot()
        self.start_animation()
        self.status_lbl.configure(text=f"Replaying: {year} {circuit}", text_color="#2CC985")

    def setup_plot(self):
        self.ax.clear()
        self.ax.axis('off')
        
        ref_driver = list(self.telemetry_data.keys())[0]
        ref_x = self.telemetry_data[ref_driver]["X"]
        ref_y = self.telemetry_data[ref_driver]["Y"]
        
        mask = ~np.isnan(ref_x)
        self.ax.plot(ref_x[mask], ref_y[mask], color='#333333', linewidth=6, alpha=0.5)
        
        # Lap Counter Text directly on Plot
        self.lap_counter_text = self.ax.text(0.02, 0.95, "Lap 1", transform=self.ax.transAxes, 
                                            color='white', fontsize=18, fontweight='bold')

        self.driver_dots = {}
        self.driver_labels = {}
        
        for driver, data in self.telemetry_data.items():
            color = data["Color"]
            dot, = self.ax.plot([], [], marker='o', color=color, markersize=8, markeredgecolor='black', markeredgewidth=1)
            self.driver_dots[driver] = dot
            
            # Label with slight transparency
            text = self.ax.text(0, 0, driver, color=color, fontsize=9, fontweight='bold', clip_on=True, 
                                bbox=dict(facecolor='black', alpha=0.5, edgecolor='none', pad=1))
            self.driver_labels[driver] = text
            
        valid_x = ref_x[mask]
        valid_y = ref_y[mask]
        self.ax.set_xlim(np.min(valid_x) - 1000, np.max(valid_x) + 1000)
        self.ax.set_ylim(np.min(valid_y) - 1000, np.max(valid_y) + 1000)
        self.ax.set_aspect('equal')
        self.canvas.draw()

    def change_speed(self, choice):
        if self.anim:
            self.start_animation(start_frame=self.current_frame)

    def start_animation(self, start_frame=0):
        if self.anim:
            try: self.anim.event_source.stop()
            except: pass
        
        total_frames = len(self.common_time)
        speed_mult = int(self.speed_var.get().replace("x", ""))
        
        # Calculate interval: 0.2s data step.
        base_interval = (0.2 * 1000) / speed_mult
        
        step_frames = 1
        if base_interval < 15: # Cap min interval to prevent UI freeze
            step_frames = int(15 / base_interval) + 1
            base_interval = 15
            
        def update(frame_idx):
            idx = start_frame + (frame_idx * step_frames)
            self.current_frame = idx
            
            if idx >= total_frames:
                try: self.anim.event_source.stop()
                except: pass
                return []

            # 1. Update Drivers & Leaderboard
            leaderboard_data = []
            
            for driver, dot in self.driver_dots.items():
                d_data = self.telemetry_data[driver]
                x, y, dist = d_data["X"][idx], d_data["Y"][idx], d_data["Distance"][idx]
                
                if np.isnan(x) or np.isnan(y):
                    dot.set_visible(False)
                    self.driver_labels[driver].set_visible(False)
                else:
                    dot.set_visible(True)
                    self.driver_labels[driver].set_visible(True)
                    dot.set_data([x], [y])
                    self.driver_labels[driver].set_position((x + 200, y + 200))
                    leaderboard_data.append((driver, dist))

            # 2. Update Lap Counter
            current_race_time = self.common_time[idx]
            lap_idx = np.searchsorted(self.lap_start_times, current_race_time, side='right') - 1
            if lap_idx >= 0 and lap_idx < len(self.lap_numbers):
                self.lap_counter_text.set_text(f"LAP {self.lap_numbers[lap_idx]} / {self.lap_numbers[-1]}")
            
            # 3. Update Leaderboard (Every 5th frame to save CPU)
            if frame_idx % 5 == 0:
                leaderboard_data.sort(key=lambda x: x[1], reverse=True) # Sort by Distance
                leader_dist = leaderboard_data[0][1] if leaderboard_data else 0
                
                for i, (row_frame, pos_lbl, drv_lbl, gap_lbl) in enumerate(self.lb_rows):
                    if i < len(leaderboard_data):
                        row_frame.pack(fill="x", pady=2) # Ensure visible
                        drv, dist = leaderboard_data[i]
                        
                        drv_lbl.configure(text=drv)
                        
                        # Style: Highlight leader
                        if i == 0:
                            gap_lbl.configure(text="LEADER", text_color="#2CC985")
                            drv_lbl.configure(text_color="#2CC985")
                        else:
                            gap = leader_dist - dist
                            gap_lbl.configure(text=f"-{gap:.0f} m", text_color="white")
                            drv_lbl.configure(text_color="white")
                    else:
                        row_frame.pack_forget() # Hide unused rows
            
            return list(self.driver_dots.values()) + list(self.driver_labels.values()) + [self.lap_counter_text]

        self.anim = FuncAnimation(self.fig, update, frames=(total_frames - start_frame) // step_frames, 
                                  interval=base_interval, blit=True, repeat=False)
        self.canvas.draw()

if __name__ == "__main__":
    app = RaceReplayerApp()
    app.mainloop()