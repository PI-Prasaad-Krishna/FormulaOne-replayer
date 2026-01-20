# F1 Race Replayer v2.5 (FastF1 + Matplotlib + CustomTkinter)
# -----------------------------------------------------------------------------
# Features:
# - MODERN GUI: CustomTkinter.
# - SMOOTH ANIMATION: 0.2s data resolution.
# - ROBUST LEADERBOARD: Uses (Laps Completed + Lap Progress) for sorting.
# - FULL RACE SUPPORT: Stitches laps for all 20 drivers.
# - PAUSE/RESUME: Functional playback control with disabled state at finish.
# - FIXED: Lap Counter timing, Leaderboard glitching & Pause Crash.
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
import customtkinter as ctk # PIP INSTALL CUSTOMTKINTER
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
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

def get_driver_color(driver_code):
    try:
        # Check if new function exists (fastf1 v3.1+)
        if hasattr(fastf1.plotting, 'get_driver_color'):
            return fastf1.plotting.get_driver_color(driver_code, session=None)
        else:
            return fastf1.plotting.driver_color(driver_code)
    except:
        return "#ffffff"

def load_race_data(year, circuit, session_type="R"):
    print(f"Loading {year} {circuit} [{session_type}] session data...")
    try:
        session = fastf1.get_session(year, circuit, session_type)
        session.load(telemetry=True, laps=True, weather=True) # EXTENDED: Load Weather
        return session
    except Exception as e:
        print(f"Error loading session: {e}")
        return None

def process_telemetry(session):
    print("Processing driver telemetry...")
    telemetry_data = {}
    drivers = pd.unique(session.laps['Driver'])
    step_size = 0.2 
    
    # Pre-calculate approximate track length from the leader's data (usually consistent)
    try:
        fastest_lap = session.laps.pick_fastest()
        track_length = fastest_lap.get_telemetry()['Distance'].max()
    except:
        track_length = 5000 # Fallback 5km

    # 1. Normalize Start Time
    # We rely on Min Telemetry Time as the absolute T=0 anchor to ensure no negative timestamps.
    # To fix the "Late Lap Counter" issue, we will calculate an offset to jump to Lap 1 immediately.
    
    min_session_time = float('inf')
    
    for driver in drivers:
        try:
            laps = session.laps.pick_drivers(driver)
            if laps.empty: continue
            tel = laps.get_telemetry()
            if tel is None or tel.empty: continue
            t_min = tel['Time'].dt.total_seconds().min()
            if t_min < min_session_time:
                min_session_time = t_min
        except: pass
        
    if min_session_time == float('inf'): return None, None, None, None, 0
    
    # We use min_session_time as the global zero
    global_start_time = min_session_time


    max_race_time = 0

    for driver in drivers:
        try:
            laps = session.laps.pick_drivers(driver)
            try:
                # Optimized approach: get full telemetry, assume continuous
                tel = laps.get_telemetry()
            except:
                continue

            if tel is None or tel.empty:
                continue
            
            # Normalize Time
            t_seconds = tel['Time'].dt.total_seconds().to_numpy() - global_start_time
            x = tel['X'].to_numpy()
            y = tel['Y'].to_numpy()
            
            # Distance logic for Leaderboard
            if 'Distance' in tel.columns:
                d = tel['Distance'].to_numpy()
            else:
                d = np.zeros_like(x) # Placeholder
            
            status = "Running"
            
            telemetry_data[driver] = {
                "Time": t_seconds,
                "X": x,
                "Y": y,
                "Distance": d, # Raw distance for interpolation
                "Color": get_driver_color(driver),
                "Team": laps.iloc[0]['Team'] if not laps.empty else "Unknown",
                "Status": status,
                "Laps": laps # Store laps for lap counting
            }
            
            if len(t_seconds) > 0:
                max_race_time = max(max_race_time, t_seconds[-1])
            
            # Robust Max Time Check using Lap Data (if available)
            try:
                if not laps.empty:
                    last_lap = laps.iloc[-1]
                    # Check for NaT (Not a Time) or NaN
                    if pd.notna(last_lap['LapStartTime']) and pd.notna(last_lap['LapTime']):
                        end_t = (last_lap['LapStartTime'] + last_lap['LapTime']).total_seconds() - global_start_time
                        if end_t > max_race_time:
                            max_race_time = end_t
            except: pass
                
        except Exception:
            pass

    if not telemetry_data:
        return None, None, None, None

    # Create Common Timeline
    common_time = np.arange(0, max_race_time, step_size)
    
    # Interpolate Data & Calculate robust metrics
    final_data = {}
    
    for driver, data in telemetry_data.items():
        # Interpolate X, Y, Distance
        
        interp_x = np.interp(common_time, data["Time"], data["X"], left=np.nan, right=data["X"][-1])
        interp_y = np.interp(common_time, data["Time"], data["Y"], left=np.nan, right=data["Y"][-1])
        interp_d = np.interp(common_time, data["Time"], data["Distance"], left=0, right=data["Distance"][-1])
        
        # Calculate End Time for this driver
        driver_end_time = data["Time"][-1]
        
        # Determine DNF status using Official Session Results
        # If driver is not classified or has a status indicating retirement
        official_status = "Unknown"
        is_dnf_overall = False
        
        try:
             # Look up driver in session.results
             # session.results is typically indexed by position or number, so search by Abbreviation
             drv_res = session.results[session.results['Abbreviation'] == driver]
             if not drv_res.empty:
                 official_status = str(drv_res.iloc[0]['Status'])
                 
                 # Define what counts as a "Finished" status
                 # 'Finished', 'Lapped', '+1 Lap', etc.
                 is_finished = (official_status.lower() in ['finished', 'lapped']) or (official_status.startswith('+'))
                 
                 is_dnf_overall = not is_finished
             else:
                 # Fallback if driver not in results (rare)
                 is_dnf_overall = driver_end_time < (max_race_time - 120)
        except Exception:
             # Fallback on error
             is_dnf_overall = driver_end_time < (max_race_time - 120)
        
        final_data[driver] = {
            "X": interp_x,
            "Y": interp_y,
            "Distance": interp_d,
            "Color": data["Color"],
            "Team": data["Team"],
            "EndTime": driver_end_time,
            "IsDNF": is_dnf_overall,
            "Laps": data["Laps"]
        }
    
    # Calculate Lap Start Times (for the Lap Counter)
    try:
        if not ref_driver:
             # Look up ref_driver again if needed (should be set above)
             max_laps = 0
             ref_driver = drivers[0]
             for d in drivers:
                try:
                    n = session.laps.pick_drivers(d)['LapNumber'].max()
                    if n > max_laps:
                        max_laps = n
                        ref_driver = d
                except: pass
            
        ref_laps = session.laps.pick_drivers(ref_driver)
        
        # LapStartTime relative to our global_start_time
        lap_start_times = ref_laps['LapStartTime'].dt.total_seconds().to_numpy() - global_start_time
        lap_numbers = ref_laps['LapNumber'].to_numpy()
        
        lap_start_times = np.nan_to_num(lap_start_times, nan=0.0)
        
        # Ensure strict sorting
        sort_idx = np.argsort(lap_start_times)
        lap_start_times = lap_start_times[sort_idx]
        lap_numbers = lap_numbers[sort_idx]
        
    except Exception as e:
        print(f"Lap counter error: {e}")
        lap_start_times = [0]
        lap_numbers = [1]
        
    # Calculate Race Start Offset (Jump to Lap 1)
    race_start_offset = 0
    try:
        # Check time of Lap 1 in our normalized timeline
        # lap_start_times[0] corresponds to the earliest lap start (should be Lap 1)
        # Note: lap_start_times is already normalized by global_start_time (min_session_time)
        
        # Determine index of "Lap 1"
        idx_l1 = np.where(lap_numbers == 1)[0]
        if len(idx_l1) > 0:
            # The start time of Lap 1 relative to T=0
            # We want to skip everything before this.
            # Clamp to 0 just in case.
            race_start_offset = max(0, lap_start_times[idx_l1[0]])
            print(f"Race Start Offset (Lap 1): {race_start_offset:.2f}s")
            
    except: pass
        

    # 2. Pre-Calculate Best Sectors & Speed Trap History per Lap
    # Structure: lap_sector_data[driver][lap_num] = {'S1': (time, color), 'S2':..., 'ST': speed}
    lap_sector_data = {d: {} for d in drivers}
    
    # Global Bests
    overall_best_s1 = float('inf')
    overall_best_s2 = float('inf')
    overall_best_s3 = float('inf')
    overall_best_st = 0
    
    # Driver Bests
    driver_best_s1 = {d: float('inf') for d in drivers}
    driver_best_s2 = {d: float('inf') for d in drivers}
    driver_best_s3 = {d: float('inf') for d in drivers}
    driver_best_st = {d: 0 for d in drivers}
    
    # Iterate through all laps sorted by time to simulate progression
    all_laps_sorted = session.laps.sort_values(by=['LapStartTime'])
    
    for idx, lap in all_laps_sorted.iterrows():
        driver = lap['Driver']
        lap_num = lap['LapNumber']
        
        # S1
        s1 = lap['Sector1Time'].total_seconds() if pd.notna(lap['Sector1Time']) else None
        s1_color = "gray"
        if s1:
            if s1 < overall_best_s1:
                overall_best_s1 = s1
                s1_color = "purple"
            elif s1 < driver_best_s1[driver]:
                driver_best_s1[driver] = s1
                s1_color = "green"
            else:
                 pass # No improvement (yellow in TV, but we keep gray/white or existing best?)
                 # Actually requirement says "Best ... times"
                 # If we show the "Best So Far", we should persist the best value.
        
        # S2
        s2 = lap['Sector2Time'].total_seconds() if pd.notna(lap['Sector2Time']) else None
        s2_color = "gray"
        if s2:
            if s2 < overall_best_s2:
                overall_best_s2 = s2
                s2_color = "purple"
            elif s2 < driver_best_s2[driver]:
                driver_best_s2[driver] = s2
                s2_color = "green"
                
        # S3
        s3 = lap['Sector3Time'].total_seconds() if pd.notna(lap['Sector3Time']) else None
        s3_color = "gray"
        if s3:
            if s3 < overall_best_s3:
                overall_best_s3 = s3
                s3_color = "purple"
            elif s3 < driver_best_s3[driver]:
                driver_best_s3[driver] = s3
                s3_color = "green"

        # Speed Trap (using SpeedST if available, or calc max)
        st = lap['SpeedST'] if 'SpeedST' in lap and pd.notna(lap['SpeedST']) else 0
        if st > overall_best_st: overall_best_st = st
        if st > driver_best_st[driver]: driver_best_st[driver] = st
        
        # Store the CURRENT BESTS for this driver at this lap
        # Requirement: "Best Sector ... times"
        # So at Lap N, we show the Best S1 they have achieved up to Lap N.
        
        cur_best_s1 = driver_best_s1[driver] if driver_best_s1[driver] != float('inf') else 0
        cur_best_s2 = driver_best_s2[driver] if driver_best_s2[driver] != float('inf') else 0
        cur_best_s3 = driver_best_s3[driver] if driver_best_s3[driver] != float('inf') else 0
        cur_best_st = driver_best_st[driver]
        
        # Determine colors for these BESTS
        # If their best == overall best, purple. Else if valid, green (since it's their best).
        # Wait, comparison with "current overall best".
        # If at Lap 10, their best is 20.0 and Overall Beest is 19.0 -> Green.
        
        c1 = "purple" if (cur_best_s1 == overall_best_s1 and cur_best_s1 > 0) else ("green" if cur_best_s1 > 0 else "gray")
        c2 = "purple" if (cur_best_s2 == overall_best_s2 and cur_best_s2 > 0) else ("green" if cur_best_s2 > 0 else "gray")
        c3 = "purple" if (cur_best_s3 == overall_best_s3 and cur_best_s3 > 0) else ("green" if cur_best_s3 > 0 else "gray")
        
        lap_sector_data[driver][lap_num] = {
            'S1': (cur_best_s1, c1),
            'S2': (cur_best_s2, c2),
            'S3': (cur_best_s3, c3),
            'ST': cur_best_st
        }

    return final_data, common_time, lap_start_times, lap_numbers, race_start_offset, session.weather_data, session.track_status, global_start_time, lap_sector_data

class RaceReplayerApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("F1 Race Replayer 2025 - Modern Edition")
        self.geometry("1400x900")
        self.protocol("WM_DELETE_WINDOW", self.on_close)

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(1, weight=1)

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
        
        self.is_paused = False
        self.is_finished = False
        self.pause_btn = ctk.CTkButton(self.control_frame, text="PAUSE", command=self.toggle_pause, fg_color="#E10600", width=80, font=("Roboto", 12, "bold"))
        self.pause_btn.pack(side="left", padx=5)
        self.pause_btn.configure(state="disabled")
        
        self.status_lbl = ctk.CTkLabel(self.control_frame, text="Ready", text_color="gray")
        self.status_lbl.pack(side="left", padx=20)

        # -- 1.1 Session Info Panel (Weather & Status) --
        self.session_info_frame = ctk.CTkFrame(self.control_frame, fg_color="transparent")
        self.session_info_frame.pack(side="right", padx=20)
        
        self.weather_lbl = ctk.CTkLabel(self.session_info_frame, text="Weather: --°C | --%", font=("Mono", 12))
        self.weather_lbl.pack(side="top", anchor="e")
        
        self.track_status_lbl = ctk.CTkLabel(self.session_info_frame, text="TRACK: GREEN", font=("Roboto", 12, "bold"), text_color="#2CC985")
        self.track_status_lbl.pack(side="bottom", anchor="e")

        # -- 2. Leaderboard (Left Side) --
        self.leaderboard_frame = ctk.CTkFrame(self, width=400, corner_radius=0) # Widened
        self.leaderboard_frame.grid(row=1, column=0, sticky="nsew", padx=0, pady=0)
        
        ctk.CTkLabel(self.leaderboard_frame, text="LEADERBOARD", font=("Roboto", 16, "bold"), text_color="#E10600").pack(pady=(15,10))
        
        header_row = ctk.CTkFrame(self.leaderboard_frame, fg_color="transparent")
        header_row.pack(fill="x", padx=10, pady=2)
        ctk.CTkLabel(header_row, text="POS", width=30, anchor="w", font=("Roboto", 10, "bold")).pack(side="left")
        ctk.CTkLabel(header_row, text="DRIVER", width=50, anchor="w", font=("Roboto", 10, "bold")).pack(side="left", padx=5)
        ctk.CTkLabel(header_row, text="TYRE", width=70, anchor="center", font=("Roboto", 10, "bold")).pack(side="left")
        ctk.CTkLabel(header_row, text="GAP", width=80, anchor="e", font=("Roboto", 10, "bold")).pack(side="left", padx=5)

        self.scroll_lb = ctk.CTkScrollableFrame(self.leaderboard_frame, fg_color="transparent")
        self.scroll_lb.pack(fill="both", expand=True, padx=5, pady=5)
        
        self.lb_rows = []
        for i in range(22):
            row_frame = ctk.CTkFrame(self.scroll_lb, height=25, fg_color=("gray85", "gray20"))
            row_frame.pack(fill="x", pady=2)
            
            pos_lbl = ctk.CTkLabel(row_frame, text=f"{i+1}", width=30, font=("Roboto", 12, "bold"))
            pos_lbl.pack(side="left", padx=5)
            
            drv_lbl = ctk.CTkLabel(row_frame, text="-", font=("Roboto", 12), width=50, anchor="w")
            drv_lbl.pack(side="left", padx=5)

            tyre_lbl = ctk.CTkLabel(row_frame, text="?", font=("Mono", 10), width=70, anchor="center")
            tyre_lbl.pack(side="left")

            gap_lbl = ctk.CTkLabel(row_frame, text="", font=("Mono", 11), width=80, anchor="e")
            gap_lbl.pack(side="left", padx=5)
            
            self.lb_rows.append((row_frame, pos_lbl, drv_lbl, tyre_lbl, gap_lbl))

        # -- 3. Map Canvas (Center/Right) --
        self.map_frame = ctk.CTkFrame(self, corner_radius=0, fg_color="#121212")
        self.map_frame.grid(row=1, column=1, columnspan=2, sticky="nsew")
        
        plt.style.use('dark_background')
        self.fig, self.ax = plt.subplots(figsize=(8, 6), facecolor='#121212')
        self.ax.set_facecolor('#121212')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.map_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.anim = None
        self.telemetry_data = {}
        self.common_time = []
        self.weather_data = None
        self.track_status_data = None
        self.global_start_time = 0
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

        self.telemetry_data, self.common_time, self.lap_start_times, self.lap_numbers, self.race_start_offset, self.weather_data, self.track_status_data, self.global_start_time, self.lap_sector_data = process_telemetry(session)
        if not self.telemetry_data:
             self.status_lbl.configure(text="No Telemetry Found", text_color="red")
             return
             
        self.setup_plot()
        self.start_animation(start_frame=int(self.race_start_offset / 0.2))
        self.pause_btn.configure(state="normal", text="PAUSE", fg_color="#E10600")
        self.is_paused = False
        self.is_finished = False
        self.status_lbl.configure(text=f"Replaying: {year} {circuit}", text_color="#2CC985")

    def setup_plot(self):
        self.ax.clear()
        self.ax.axis('off')
        
        ref_driver = list(self.telemetry_data.keys())[0]
        ref_x = self.telemetry_data[ref_driver]["X"]
        ref_y = self.telemetry_data[ref_driver]["Y"]
        
        mask = ~np.isnan(ref_x)
        self.ax.plot(ref_x[mask], ref_y[mask], color='#333333', linewidth=6, alpha=0.5)
        
        self.lap_counter_text = self.ax.text(0.02, 0.95, "Lap 1", transform=self.ax.transAxes, 
                                            color='white', fontsize=18, fontweight='bold')

        self.driver_dots = {}
        self.driver_labels = {}
        
        for driver, data in self.telemetry_data.items():
            color = data["Color"]
            dot, = self.ax.plot([], [], marker='o', color=color, markersize=8, markeredgecolor='black', markeredgewidth=1)
            self.driver_dots[driver] = dot
            
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
        if self.anim and not self.is_paused and not self.is_finished:
            self.start_animation(start_frame=self.current_frame)

    def toggle_pause(self):
        # FIX: Check if anim exists before trying to access event_source
        if self.anim is None or self.is_finished:
            return
        
        if self.is_paused:
            try:
                self.anim.event_source.start()
                self.pause_btn.configure(text="PAUSE", fg_color="#E10600")
                self.is_paused = False
            except AttributeError:
                pass
        else:
            try:
                self.anim.event_source.stop()
                self.pause_btn.configure(text="RESUME", fg_color="#2CC985")
                self.is_paused = True
            except AttributeError:
                pass

    def start_animation(self, start_frame=0):
        if self.anim:
            try: self.anim.event_source.stop()
            except: pass
        
        total_frames = len(self.common_time)
        speed_mult = int(self.speed_var.get().replace("x", ""))
        base_interval = (0.2 * 1000) / speed_mult
        
        step_frames = 1
        if base_interval < 15:
            step_frames = int(15 / base_interval) + 1
            base_interval = 15
            
        def update(frame_idx):
            idx = start_frame + (frame_idx * step_frames)
            self.current_frame = idx
            
            if idx >= total_frames:
                try: self.anim.event_source.stop()
                except: pass
                self.status_lbl.configure(text="Race Finished", text_color="white")
                self.pause_btn.configure(state="disabled", text="FINISHED", fg_color="gray")
                self.is_finished = True
                return []

            leaderboard_data = []
            current_race_time = self.common_time[idx]

            for driver, dot in self.driver_dots.items():
                d_data = self.telemetry_data[driver]
                
                # Check for Albon issue (or any driver with missing data)
                if idx >= len(d_data["X"]):
                    continue # Skip if index out of bounds for this driver

                x, y, dist = d_data["X"][idx], d_data["Y"][idx], d_data["Distance"][idx]
                end_time = d_data.get("EndTime", 999999)
                
                # Determine DNF status at current frame
                # If current time > end_time + buffer, they are stopped.
                is_dnf = (d_data["IsDNF"] and current_race_time > (end_time + 10))

                if np.isnan(x) or np.isnan(y):
                    dot.set_visible(False)
                    self.driver_labels[driver].set_visible(False)
                else:
                    dot.set_visible(True)
                    self.driver_labels[driver].set_visible(True)
                    dot.set_data([x], [y])
                    self.driver_labels[driver].set_position((x + 200, y + 200))
                    
                    if is_dnf:
                        dot.set_color("gray")
                        self.driver_labels[driver].set_color("gray")
                        self.driver_labels[driver].set_text(f"{driver} (DNF)")
                    else:
                        pass
                    
                    leaderboard_data.append((driver, dist, is_dnf))

            # Update Lap Counter - Fixed Logic
            lap_idx = np.searchsorted(self.lap_start_times, current_race_time, side='right') - 1
            
            if lap_idx >= 0:
                if lap_idx < len(self.lap_numbers):
                    current_lap_num = self.lap_numbers[lap_idx]
                    total_laps = self.lap_numbers[-1]
                    self.lap_counter_text.set_text(f"LAP {current_lap_num} / {total_laps}")
                else:
                    self.lap_counter_text.set_text("FINISH")
            else:
                self.lap_counter_text.set_text("GRID")
            
            # Update Session Info (Weather & Status) - Every 1 second (approx 5 frames)
            if frame_idx % 5 == 0:
                try:
                    current_session_time = pd.Timedelta(seconds=(current_race_time + self.global_start_time))
                    
                    # Weather
                    if self.weather_data is not None:
                        # Find closest weather row before current time
                        w_row = self.weather_data[self.weather_data['Time'] <= current_session_time].iloc[-1]
                        air_temp = w_row['AirTemp']
                        humidity = w_row['Humidity']
                        self.weather_lbl.configure(text=f"Air: {air_temp}°C | Hum: {humidity}%")
                    
                    # Track Status
                    if self.track_status_data is not None:
                        # Status is often '1' (Green), '2' (Yellow), '4' (SC), '5' (Red), '6' (VSC), '7' (VSC ending)
                        ts_row = self.track_status_data[self.track_status_data['Time'] <= current_session_time].iloc[-1]
                        status_code = ts_row['Status']
                        
                        status_text = "GREEN"
                        status_color = "#2CC985"
                        
                        if status_code == '1': 
                            status_text = "GREEN"
                            status_color = "#2CC985"
                        elif status_code == '2':
                            status_text = "YELLOW"
                            status_color = "#E10600" # Reddish for visibility (or yellow)
                        elif status_code == '4':
                            status_text = "SC"
                            status_color = "orange"
                        elif status_code in ['6', '7']:
                            status_text = "VSC"
                            status_color = "orange"
                        elif status_code == '5':
                            status_text = "RED FLAG"
                            status_color = "red"
                            
                        self.track_status_lbl.configure(text=f"TRACK: {status_text}", text_color=status_color)
                        
                except: pass

            # Update Leaderboard (Every 5th frame)
            if frame_idx % 5 == 0:
                # SORTING LOGIC: Total Distance
                leaderboard_data.sort(key=lambda x: x[1], reverse=True)
                
                active_leaders = [d for d in leaderboard_data if not d[2]]
                leader_dist = active_leaders[0][1] if active_leaders else 0
                
                for i, (row_frame, pos_lbl, drv_lbl, tyre_lbl, gap_lbl) in enumerate(self.lb_rows):
                    if i < len(leaderboard_data):
                        row_frame.pack(fill="x", pady=2)
                        drv, dist, is_dnf = leaderboard_data[i]
                        
                        drv_lbl.configure(text=drv)
                        
                        # Tyre Logic
                        tyre_text = "?"
                        tyre_color = "white"
                        try:
                            # Find current lap for driver
                            d_laps = self.telemetry_data[drv]['Laps']
                            
                            compound = "Unknown"
                            # Filter laps that have started
                            current_session_time = pd.Timedelta(seconds=(current_race_time + self.global_start_time))
                            
                            started_laps = d_laps[d_laps['LapStartTime'] <= current_session_time]
                            if not started_laps.empty:
                                current_lap = started_laps.iloc[-1]
                                compound = str(current_lap['Compound']).upper()
                            else:
                                # Fallback to first lap if race hasn't started for them
                                if not d_laps.empty:
                                    current_lap = d_laps.iloc[0]
                                    compound = str(current_lap['Compound']).upper()

                            if 'SOFT' in compound:
                                tyre_text = "Soft"
                                tyre_color = "#FF3333" # Red
                            elif 'MEDIUM' in compound:
                                tyre_text = "Medium"
                                tyre_color = "#FFE000" # Yellow
                            elif 'HARD' in compound:
                                tyre_text = "Hard"
                                tyre_color = "white"
                            elif 'INTER' in compound:
                                tyre_text = "Inter"
                                tyre_color = "#39B54A" # Green
                            elif 'WET' in compound:
                                tyre_text = "Wet"
                                tyre_color = "#00AEEF" # Cyan
                            else:
                                # Show first 3 chars if unknown/NAN
                                tyre_text = compound[:3] if compound and compound != "NAN" else "?"
                                tyre_color = "gray"
                            
                        except Exception as e: 
                            print(f"Tyre Error {drv}: {e}")
                            print(f"Columns: {d_laps.columns}")
                            tyre_text = "Err"
                        
                        tyre_lbl.configure(text=tyre_text, text_color=tyre_color)

                        if is_dnf:
                             gap_lbl.configure(text="DNF", text_color="red")
                             drv_lbl.configure(text_color="gray")
                             pos_lbl.configure(text_color="gray")
                        else:
                            if i == 0:
                                gap_lbl.configure(text="LEADER", text_color="#2CC985")
                                drv_lbl.configure(text_color="#2CC985")
                                pos_lbl.configure(text_color="#2CC985")
                            else:
                                gap = leader_dist - dist
                                gap_lbl.configure(text=f"-{gap:.0f} m", text_color="white")
                                drv_lbl.configure(text_color="white")
                                pos_lbl.configure(text_color="white")
                    else:
                        row_frame.pack_forget()
            
            return list(self.driver_dots.values()) + list(self.driver_labels.values()) + [self.lap_counter_text]

        self.anim = FuncAnimation(self.fig, update, frames=(total_frames - start_frame) // step_frames, 
                                  interval=base_interval, blit=True, repeat=False)
        self.canvas.draw()

if __name__ == "__main__":
    app = RaceReplayerApp()
    app.mainloop()