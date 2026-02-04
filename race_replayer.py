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
import threading

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

def format_time(seconds):
    if seconds is None or seconds < 0:
        return "00:00.000"
    minutes = int(seconds // 60)
    secs = seconds % 60
    return f"{minutes:02}:{secs:06.3f}"

def get_driver_color(driver_code, session=None):
    try:
        # Check if new function exists (fastf1 v3.1+)
        if hasattr(fastf1.plotting, 'get_driver_color'):
            return fastf1.plotting.get_driver_color(driver_code, session=session)
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
                "X": x,
                "Y": y,
                "Z": tel['Z'].to_numpy() if 'Z' in tel.columns else np.zeros_like(x), # Extract Z (Elevation)
                "Distance": d, # Raw distance for interpolation
                "Speed": tel['Speed'].to_numpy() if 'Speed' in tel.columns else np.zeros_like(x),
                "RPM": tel['RPM'].to_numpy() if 'RPM' in tel.columns else np.zeros_like(x),
                "Gear": tel['nGear'].to_numpy() if 'nGear' in tel.columns else np.zeros_like(x),
                "Throttle": tel['Throttle'].to_numpy() if 'Throttle' in tel.columns else np.zeros_like(x),
                "Brake": tel['Brake'].to_numpy() if 'Brake' in tel.columns else np.zeros_like(x),
                "Color": get_driver_color(driver, session=session),
                "Team": laps.iloc[0]['Team'] if not laps.empty else "Unknown",
                "Status": status,
                "Laps": laps, # Store laps for lap counting
                "PitIntervals": [] # Initialize
            }
            
            # Pre-calc normalized lap start times for this driver
            try:
                # Avoid modifying the original session.laps slice in place if it affects others, 
                # but pandas usually copies on filter. Safe to set here.
                # Use .loc to avoid SettingWithCopyWarning
                laps = laps.copy()
                laps['NormLapStartTime'] = laps['LapStartTime'].dt.total_seconds() - global_start_time
                
                # FIX: Fallback for Out Laps where LapStartTime might be NaT
                # Use PitOutTime if available
                if 'PitOutTime' in laps.columns:
                     pit_out_norm = laps['PitOutTime'].dt.total_seconds() - global_start_time
                     laps['NormLapStartTime'] = laps['NormLapStartTime'].fillna(pit_out_norm)
                
                # Final Fallback: Forward Fill (Propagate previous lap start)
                # This ensures we don't drop the lap entirely, even if timing is imprecise
                laps['NormLapStartTime'] = laps['NormLapStartTime'].ffill()
                
                laps = laps.sort_values(by='LapNumber') # Enforce sorting
                telemetry_data[driver]["Laps"] = laps
            except Exception as e:
                print(f"Error pre-calculating lap times for {driver}: {e}")
            
            # Calculate Pit Intervals (Robust)
            try:
                # Get all valid pit times
                pit_ins = laps['PitInTime']
                pit_outs = laps['PitOutTime']
                
                # Convert to seconds (handle NaT safely)
                raw_ins = pit_ins.dropna().dt.total_seconds().values
                raw_outs = pit_outs.dropna().dt.total_seconds().values
                
                p_intervals = []
                BUFFER = 3.0 # Increased to 3.0s to ensure coverage, reliance on timing only
                
                # 1. Start from Pit Lane (Out without preceding In)
                if len(raw_outs) > 0:
                    first_out = raw_outs[0]
                    # Only calculate if Out < In (or no In) AND it happens early (first 15 mins)
                    # Otherwise it's just a normal stop where we missed the In time (rare) or a very long first stint
                    is_early_start = (first_out - global_start_time) < 900 
                    
                    if (len(raw_ins) == 0 or first_out < raw_ins[0]) and is_early_start:
                        # Driver started in pits, so from T=0 to PitOut
                        p_intervals.append( (-10, first_out - global_start_time + BUFFER) )
                
                # 2. Standard Pit Stops
                for t_in in raw_ins:
                     start_t = t_in - global_start_time - BUFFER
                     
                     # Find corresponding out
                     valid_outs = raw_outs[raw_outs > t_in]
                     if len(valid_outs) > 0:
                         end_t = valid_outs[0] - global_start_time + BUFFER
                     else:
                         # Retirement in pits
                         end_t = start_t + 120
                         
                     p_intervals.append( (start_t, end_t) )
                
                telemetry_data[driver]["PitIntervals"] = p_intervals
            except Exception as e:
                print(f"Pit interval calc error {driver}: {e}")
                telemetry_data[driver]["PitIntervals"] = []
            
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
        return None, None, None, None, None, None, None, None, None, None

    # Create Common Timeline
    common_time = np.arange(0, max_race_time, step_size)
    
    # Interpolate Data & Calculate robust metrics
    final_data = {}
    
    for driver, data in telemetry_data.items():
        # Interpolate X, Y, Distance
        
        interp_x = np.interp(common_time, data["Time"], data["X"], left=np.nan, right=data["X"][-1])
        interp_y = np.interp(common_time, data["Time"], data["Y"], left=np.nan, right=data["Y"][-1])
        interp_z = np.interp(common_time, data["Time"], data["Z"], left=np.nan, right=data["Z"][-1])
        interp_d = np.interp(common_time, data["Time"], data["Distance"], left=0, right=data["Distance"][-1])
        interp_s = np.interp(common_time, data["Time"], data["Speed"], left=0, right=0) 
        interp_rpm = np.interp(common_time, data["Time"], data["RPM"], left=0, right=0)
        interp_gear = np.interp(common_time, data["Time"], data["Gear"], left=0, right=0)
        interp_thr = np.interp(common_time, data["Time"], data["Throttle"], left=0, right=0)
        interp_brk = np.interp(common_time, data["Time"], data["Brake"], left=0, right=0)

        
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
            "Z": interp_z,
            "Distance": interp_d,
            "Speed": interp_s,
            "RPM": interp_rpm,
            "Gear": interp_gear,
            "Throttle": interp_thr,
            "Brake": interp_brk,
            "Color": data["Color"],
            "Team": data["Team"],
            "EndTime": driver_end_time,
            "IsDNF": is_dnf_overall,
            "Laps": data["Laps"],
            "PitIntervals": data.get("PitIntervals", [])
        }
    
    # 3. Build Pit Lane Path (Geometry) for Proximity Check
    # Collect X,Y points from all drivers where they are in a known Pit Interval
    pit_lane_points_x = []
    pit_lane_points_y = []
    
    BUFFER = 1.0 # Reduced from 8.0 to 1.0 match above
    
    for driver, data in final_data.items():
        intervals = data["PitIntervals"]
        for (start_t, end_t) in intervals:
            # STRICT MODE + SPEED FILTER
            # Exclude buffer from time range
            strict_start = start_t + BUFFER
            strict_end = end_t - BUFFER
            
            if strict_end <= strict_start: continue
            
            # Find indices in common_time
            mask = (common_time >= strict_start) & (common_time <= strict_end)
            if np.any(mask):
                 # Get X, Y, Speed
                 px = data["X"][mask]
                 py = data["Y"][mask]
                 ps = data["Speed"][mask]
                 
                 # ONLY add points where Speed < 90 km/h (actual pit lane travel)
                 # This filters out high-speed entry/exit ramps that might overlap with track
                 speed_mask = ps < 90
                 
                 if np.any(speed_mask):
                     pit_lane_points_x.extend(px[speed_mask])
                     pit_lane_points_y.extend(py[speed_mask])
    
    pit_lane_path = None
    if pit_lane_points_x:
        # Downsample
        pts = np.column_stack((pit_lane_points_x, pit_lane_points_y))
        if len(pts) > 0:
            pit_lane_path = pts[::5] # Store as numpy array
            
    # Calculate Lap Start Times (for the Lap Counter)
    try:
        # Identify reference driver (most laps) to validat ref_driver variable
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
    # Calculate Race Start Offset via MOVEMENT DETECTION
    # The Lap 1 Start Time metadata seems unreliable (pointing to mid-race or restart).
    # We will find the first time any driver exceeds 10 km/h (Formation Start).
    race_start_offset = 0
    try:
        min_move_time = float('inf')
        
        # Check a few top drivers to avoid outliers
        check_drivers = drivers[:5] if len(drivers) > 5 else drivers
        
        for d in check_drivers:
            if d in telemetry_data:
                d_obj = telemetry_data[d]
                # Find first index where speed > 10
                move_idx = np.argmax(d_obj["Speed"] > 10)
                if d_obj["Speed"][move_idx] > 10: # verify it actually found one
                    t_move = d_obj["Time"][move_idx]
                    if t_move < min_move_time:
                        min_move_time = t_move
                        
        if min_move_time != float('inf'):
            # Start 2 mins before movement (Formation Lap buffer)
            race_start_offset = max(0, min_move_time - 120)
            print(f"Detected First Movement: {min_move_time:.2f}s")
            print(f"Race Start Offset (Movement - 120s): {race_start_offset:.2f}s")
        else:
            print("No movement detected, starting at 0.")

    except Exception as e:
        print(f"Movement Calc Error: {e}")
        

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

    return final_data, common_time, lap_start_times, lap_numbers, race_start_offset, session.weather_data, session.track_status, global_start_time, lap_sector_data, pit_lane_path

class RaceReplayerApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("F1 Race Replayer 2025 - Modern Edition")
        self.geometry("1400x900")
        self.protocol("WM_DELETE_WINDOW", self.on_close)
        
        # State for Selection
        self.selected_driver = None
        self.leaderboard_data_ref = [] # To store current frame's leaderboard order
        
        # Animation State
        self.anim = None



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
        
        # 3D Toggle Switch
        self.is_3d_mode = ctk.BooleanVar(value=False)
        self.switch_3d = ctk.CTkSwitch(self.control_frame, text="3D MAP", variable=self.is_3d_mode, command=self.toggle_3d_mode, progress_color="#E10600")
        self.switch_3d.pack(side="left", padx=20)

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
            
            # Make interactive
            # We use a closure to capture 'i' correctly
            def make_handler(index):
                return lambda e: self.on_row_click(index)
            
            handler = make_handler(i)
            
            row_frame.bind("<Button-1>", handler)
            pos_lbl.bind("<Button-1>", handler)
            drv_lbl.bind("<Button-1>", handler)
            tyre_lbl.bind("<Button-1>", handler)
            gap_lbl.bind("<Button-1>", handler)
            
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
        self.hud_elements = {} # Store HUD artists
        self.current_frame = 0
        
        # Resize Handling
        self.resize_timer = None
        self.was_playing = False
        self.map_frame.bind("<Configure>", self.on_resize)

        # Bind scroll event for 3D zoom
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

    def on_row_click(self, row_idx):
        """ Handle click on leaderboard row """
        if not self.leaderboard_data_ref or row_idx >= len(self.leaderboard_data_ref):
            return
            
        clicked_driver = self.leaderboard_data_ref[row_idx][0]
        
        if self.selected_driver == clicked_driver:
            # Deselect
            self.selected_driver = None
            print(f"[Selection] Deselected {clicked_driver}. Returning to Leader.")
        else:
            # Select
            self.selected_driver = clicked_driver
            print(f"[Selection] Selected {clicked_driver}")
            
        # Force update will happen next frame, or we could force it here if paused.
        # If paused, we might want to manually refresh HUD:
        if self.is_paused:
             # Basic refresh if possible, but data might be stale. 
             # Simplest is just to wait for next frame or user resume.
             pass

    def on_scroll(self, event):
        """ Zoom in/out with scroll wheel in 3D mode by adjusting axis limits """
        if not self.is_3d_mode.get():
            return
            
        if event.inaxes != self.ax:
            return

        # Zoom Factor logic
        # < 1.0 = Zoom In (Smaller View)
        # > 1.0 = Zoom Out (Larger View)
        zoom_step = 0.1
        
        if event.button == 'up':
            # Zoom In
            self.zoom_level = max(0.1, self.zoom_level - zoom_step)
        elif event.button == 'down':
            # Zoom Out
            self.zoom_level = min(5.0, self.zoom_level + zoom_step)
        else:
            return
            
        self.apply_zoom()
        self.canvas.draw_idle()

    def apply_zoom(self):
        """ Apply current zoom level to axis limits """
        if not hasattr(self, 'plot_bounds'):
            return
            
        (min_x, max_x, min_y, max_y, min_z, max_z) = self.plot_bounds
        
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        
        span_x = (max_x - min_x)
        span_y = (max_y - min_y)
        
        # Apply Zoom
        new_span_x = span_x * self.zoom_level
        new_span_y = span_y * self.zoom_level
        
        self.ax.set_xlim(center_x - new_span_x/2, center_x + new_span_x/2)
        self.ax.set_ylim(center_y - new_span_y/2, center_y + new_span_y/2)
        
        # Optionally zoom Z or keep it fixed? 
        # Usually for track maps we just zoom X/Y. 
        # But if we zoom in a lot, Z might need adjustment if we are looking at elevation changes.
        # Let's keep Z fixed for now or maybe zoom it less aggressively.
        # self.ax.set_zlim(...) 

    def on_close(self):
        if getattr(self, 'anim', None):
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

        # Disable Load Button to prevent multiple clicks
        self.load_btn.configure(state="disabled", text="LOADING...")
        self.status_lbl.configure(text="Loading Data... (Please Wait)", text_color="orange")
        
        # Stop existing animation
        if self.anim: 
             try: self.anim.event_source.stop()
             except: pass
        self.update() # Force UI update immediately

        # Start Thread
        loader_thread = threading.Thread(target=self.perform_loading, args=(year, circuit), daemon=True)
        loader_thread.start()

    def perform_loading(self, year, circuit):
        """ Run in background thread """
        try:
            session = load_race_data(year, circuit, "R")
            if not session:
                self.after(0, self.finish_loading, (False, "Could not load session data."))
                return

            results = process_telemetry(session)
            # results might be (None, ...) if telemetry failed, check first element
            if not results or results[0] is None:
                 self.after(0, self.finish_loading, (False, "No Telemetry Found"))
            else:
                 self.after(0, self.finish_loading, (True, results))
                 
        except Exception as e:
            self.after(0, self.finish_loading, (False, f"Error: {str(e)}"))

    def finish_loading(self, result):
        """ Run in main thread """
        success, data = result
        
        # Re-enable button
        self.load_btn.configure(state="normal", text="LOAD RACE")

        if not success:
            error_msg = data
            self.status_lbl.configure(text=f"Load Failed: {error_msg}", text_color="red")
            messagebox.showerror("Error", error_msg)
            return
        
        # Unpack Data
        # (telemetry_data, common_time, lap_start_times, lap_numbers, race_start_offset, 
        #  weather_data, track_status_data, global_start_time, lap_sector_data, pit_lane_path)
        (self.telemetry_data, self.common_time, self.lap_start_times, self.lap_numbers, 
         self.race_start_offset, self.weather_data, self.track_status_data, 
         self.global_start_time, self.lap_sector_data, self.pit_lane_path) = data
             
        self.setup_plot()
        self.start_animation(start_frame=int(self.race_start_offset / 0.2))
        self.pause_btn.configure(state="normal", text="PAUSE", fg_color="#E10600")
        self.is_paused = False
        self.is_finished = False
        
        year = self.year_entry.get()
        circuit = self.circuit_entry.get()
        self.status_lbl.configure(text=f"Replaying: {year} {circuit}", text_color="#2CC985")

                                    
    def toggle_3d_mode(self):
        # Restart plot with new mode if data exists
        if self.telemetry_data:
            self.setup_plot()
            if not self.is_paused and not self.is_finished:
                 self.start_animation(start_frame=self.current_frame)
            else:
                 # Just redraw current frame static
                 pass

    def setup_plot(self):
        # safely clear
        self.fig.clf() # Clear figure to reset axes (important for 2D vs 3D)
        
        # Check Mode
        use_3d = self.is_3d_mode.get()
        
        # --- 1. SPLIT FIGURES (Map vs HUD) ---
        # Map: Left 75%
        # HUD: Right 25%
        
        if use_3d:
            self.ax = self.fig.add_axes([0.0, 0.0, 0.75, 1.0], projection='3d')
            self.ax.set_facecolor('#121212')
            # Pane colors
            self.ax.xaxis.set_pane_color((0.1, 0.1, 0.1, 1.0))
            self.ax.yaxis.set_pane_color((0.1, 0.1, 0.1, 1.0))
            self.ax.zaxis.set_pane_color((0.1, 0.1, 0.1, 1.0))
            self.ax.grid(False)
        else:
            self.ax = self.fig.add_axes([0.0, 0.0, 0.75, 1.0])
            self.ax.axis('off')
            self.ax.set_facecolor('#121212')
            
        # HUD Axes (Always 2D)
        self.ax_hud = self.fig.add_axes([0.75, 0.0, 0.25, 1.0])
        self.ax_hud.set_facecolor('#181818') # Slightly lighter background for sidebar
        self.ax_hud.axis('off')
        
        # Add visual separator line
        self.ax_hud.axvline(0, color='#333333', linewidth=2)
        
        ref_driver = list(self.telemetry_data.keys())[0]
        ref_x = self.telemetry_data[ref_driver]["X"]
        ref_y = self.telemetry_data[ref_driver]["Y"]
        ref_z = self.telemetry_data[ref_driver]["Z"]
        
        mask = ~np.isnan(ref_x)
        
        if use_3d:
             self.ax.plot(ref_x[mask], ref_y[mask], ref_z[mask], color='#333333', linewidth=4, alpha=0.5)
        else:
             self.ax.plot(ref_x[mask], ref_y[mask], color='#333333', linewidth=6, alpha=0.5)
        
        # Lap Counter (Top Left of Map)
        if use_3d:
            self.lap_counter_text = self.ax.text2D(0.02, 0.95, "00:00.000", transform=self.ax.transAxes, 
                                                color='white', fontsize=18, fontweight='bold')
        else:
            self.lap_counter_text = self.ax.text(0.02, 0.95, "00:00.000", transform=self.ax.transAxes, 
                                                color='white', fontsize=18, fontweight='bold')

        self.driver_dots = {}
        self.driver_labels = {}
        self.hud_elements = {} 
        
        for driver, data in self.telemetry_data.items():
            color = data["Color"]
            if use_3d:
                 # 3D Plot
                 dot, = self.ax.plot([], [], [], marker='o', color=color, markersize=6, markeredgecolor='black', markeredgewidth=1)
            else:
                 dot, = self.ax.plot([], [], marker='o', color=color, markersize=8, markeredgecolor='black', markeredgewidth=1)
            
            self.driver_dots[driver] = dot
            
            # Text in 3D
            if use_3d:
                 text = self.ax.text(0, 0, 0, driver, color=color, fontsize=8, fontweight='bold')
            else:
                 text = self.ax.text(0, 0, driver, color=color, fontsize=9, fontweight='bold', clip_on=True, 
                                bbox=dict(facecolor='black', alpha=0.5, edgecolor='none', pad=1))
            self.driver_labels[driver] = text
        
        # --- HUD LAYOUT (Sidebar) ---
        
        # Helper for HUD text (uses ax_hud)
        def add_hud_text(x, y, s, **kwargs):
             return self.ax_hud.text(x, y, s, transform=self.ax_hud.transAxes, **kwargs)

        # 1. HEADER (Top)
        add_hud_text(0.5, 0.92, "TELEMETRY", color='gray', fontsize=10, ha='center', weight='bold')
        self.hud_driver = add_hud_text(0.5, 0.88, "LEADER", color='white', 
                                  fontsize=14, fontweight='bold', ha='center',
                                  bbox=dict(facecolor='#222222', alpha=0.8, pad=4, edgecolor='none'))

        # 2. SPEED (Upper Middle)
        self.hud_speed = add_hud_text(0.5, 0.70, "0", color='white', fontsize=40, fontweight='bold', ha='center')
        self.hud_unit = add_hud_text(0.5, 0.64, "km/h", color='gray', fontsize=12, ha='center')
        
        # 3. GEAR & RPM (Middle)
        # Gear Circle background
        # self.ax_hud.add_patch(plt.Circle((0.5, 0.52), 0.08, transform=self.ax_hud.transAxes, color='#222222'))
        self.hud_gear = add_hud_text(0.5, 0.52, "N", color='#00ff00', fontsize=30, fontweight='bold', ha='center', va='center')
        self.hud_rpm = add_hud_text(0.5, 0.45, "0 RPM", color='lightgray', fontsize=12, ha='center')
        
        # 4. BARS (Bottom)
        # Vertical arrangement for sidebar? Or Horizontal?
        # Let's do Horizontal bars but stacked.
        
        # Throttle
        add_hud_text(0.1, 0.30, "THR", color='gray', fontsize=10, ha='left')
        self.hud_thr_text = add_hud_text(0.9, 0.30, "0%", color='#2CC985', fontsize=10, ha='right')
        
        self.ax_hud.add_patch(plt.Rectangle((0.10, 0.27), 0.80, 0.02, transform=self.ax_hud.transAxes, color='#333333'))
        self.hud_thr_bar = plt.Rectangle((0.10, 0.27), 0.0, 0.02, transform=self.ax_hud.transAxes, color='#2CC985')
        self.ax_hud.add_patch(self.hud_thr_bar)
        
        # Brake
        add_hud_text(0.1, 0.20, "BRK", color='gray', fontsize=10, ha='left')
        self.hud_brk_text = add_hud_text(0.9, 0.20, "0%", color='#E10600', fontsize=10, ha='right')
        
        self.ax_hud.add_patch(plt.Rectangle((0.10, 0.17), 0.80, 0.02, transform=self.ax_hud.transAxes, color='#333333'))
        self.hud_brk_bar = plt.Rectangle((0.10, 0.17), 0.0, 0.02, transform=self.ax_hud.transAxes, color='#E10600')
        self.ax_hud.add_patch(self.hud_brk_bar)

        valid_x = ref_x[mask]
        valid_y = ref_y[mask]
        valid_z = ref_z[mask]
        
        # Store bounds for zooming
        min_x, max_x = np.min(valid_x) - 500, np.max(valid_x) + 500
        min_y, max_y = np.min(valid_y) - 500, np.max(valid_y) + 500
        min_z, max_z = np.min(valid_z) - 50, np.max(valid_z) + 50
        
        self.plot_bounds = (min_x, max_x, min_y, max_y, min_z, max_z)
        self.zoom_level = 1.0
        
        if use_3d:
             self.ax.set_zlim(min_z, max_z)
             # Set view
             self.ax.view_init(elev=30, azim=-60)
             self.ax.set_box_aspect((1, 1, 0.2)) # Flatter Z
        else:
             self.ax.set_aspect('equal')
             
        self.apply_zoom()
        self.canvas.draw()
        
    def on_resize(self, event):
        """ Debounce resize events to prevent crash with blit=True """
        # Only care if map frame size changed
        if self.resize_timer:
            self.after_cancel(self.resize_timer)
            
        if self.anim and self.anim.event_source:
             try:
                 self.anim.event_source.stop()
                 if not self.is_paused and not self.is_finished:
                     self.was_playing = True
                     self.is_paused = True # Mark as paused temporarily
             except: pass
        
        # Schedule finish_resize
        self.resize_timer = self.after(500, self.finish_resize)

    def finish_resize(self):
        """ Restart animation after resize settles """
        self.canvas.draw() # Redraw to fix resolution
        
        if self.was_playing:
            self.is_paused = False
            self.was_playing = False
            self.start_animation(self.current_frame)
        elif self.telemetry_data:
             # Just redraw static if we were paused
             # self.setup_plot() # Might be too heavy, but draw() above captures it
             pass

    def change_speed(self, choice):
        if self.anim and not self.is_paused and not self.is_finished:
            self.start_animation(start_frame=self.current_frame)

    def toggle_pause(self):
        # FIX: Check if anim exists before trying to access event_source
        if self.anim is None or self.is_finished:
            return
        
        if self.is_paused:
            try:
                # RESUME: Call start_animation to pick up new speed settings
                self.start_animation(start_frame=self.current_frame)
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
                # FIX: Instead of continue (which skips update), clamp to last frame
                # This ensures the dot stays at the final position if data runs out slightly early
                d_len = len(d_data["X"])
                if d_len == 0: continue
                
                safe_idx = min(idx, d_len - 1)

                x, y, z, dist = d_data["X"][safe_idx], d_data["Y"][safe_idx], d_data["Z"][safe_idx], d_data["Distance"][safe_idx]
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
                    
                    if self.is_3d_mode.get():
                        dot.set_data([x], [y])
                        dot.set_3d_properties([z])
                        self.driver_labels[driver].set_position((x, y))
                        self.driver_labels[driver].set_3d_properties(z)
                    else:
                        dot.set_data([x], [y])
                        self.driver_labels[driver].set_position((x + 200, y + 200))
                    
                    # 1. Check Pit Status (Timing + Proximity Fallback)
                    is_pitting = False
                    
                    # A. Timing Interval Check
                    if "PitIntervals" in d_data:
                        for (p_start, p_end) in d_data["PitIntervals"]:
                            if p_start <= current_race_time <= p_end:
                                is_pitting = True
                                break
                    
                    # B. Proximity Fallback Check (if not already detected)
                    # B. Proximity Fallback Check REMOVED
                    # We rely solely on official start/end timing intervals to avoid false positives on track corners.
                    pass

                    dot.set_data([x], [y])
                    self.driver_labels[driver].set_position((x + 200, y + 200))
                    
                    if is_dnf:
                         label_text = f"{driver} (DNF)"
                    else:
                         label_text = driver

                    if is_dnf:
                         dot.set_color("gray")
                         self.driver_labels[driver].set_color("gray")
                         self.driver_labels[driver].set_text(label_text)
                    else:
                        # Restore color
                        dot.set_color(d_data["Color"])
                        self.driver_labels[driver].set_color(d_data["Color"])
                        self.driver_labels[driver].set_text(label_text)
                    
                    leaderboard_data.append((driver, dist, is_dnf))

            # Update Lap Counter - Dynamic FastF1 Logic
            try:
                # 1. Identify Leader (Driver with max Distance)
                # Filter out DNFs from potential leaders to avoid stuck counters if leader retires
                active_runners = [d for d in leaderboard_data if not d[2]] 
                if not active_runners:
                     active_runners = leaderboard_data # Fallback if everyone DNF (unlikely)

                if active_runners:
                    # Sort by distance descending
                    active_runners.sort(key=lambda x: x[1], reverse=True)
                    leader_driver = active_runners[0][0]
                    
                    # Define HUD Driver immediately available for use
                    hud_driver_name = self.selected_driver if self.selected_driver else leader_driver
                    
                    # 2. Get Leader's Lap Data
                    if leader_driver in self.telemetry_data:
                        laps_df = self.telemetry_data[leader_driver]["Laps"]
                        
                        # 3. Find Current Lap based on Time
                        # current_race_time is normalized (0 start). 
                        # We need to compare with LapStartTime normalized.
                        
                        # We can't vectorise easily here without pre-calc, but iterating 60 rows is fine.
                        # We look for the last lap where LapStartTime <= current_race_time
                        
                        current_lap_num = 1
                        total_laps = 0
                        
                        if not laps_df.empty:
                            total_laps = int(laps_df['LapNumber'].max())
                            
                            # Filter laps that have started
                            # We must handle NaT or invalid times carefully
                            valid_laps = laps_df[pd.notna(laps_df['LapStartTime'])].copy()
                            
                            if not valid_laps.empty:
                                # Calculate normalized start times for comparison
                                # NOTE: This could be optimized by pre-calculating in process_telemetry, 
                                # but it's fast enough here.
                                valid_laps['NormStartTime'] = valid_laps['LapStartTime'].dt.total_seconds() - self.global_start_time
                                
                                # Find laps started before now
                                started_laps = valid_laps[valid_laps['NormStartTime'] <= current_race_time]
                                
                                if not started_laps.empty:
                                    current_lap_num = int(started_laps.iloc[-1]['LapNumber'])
                                else:
                                    current_lap_num = 0

                        # Update Timer instead of Lap Counter
                        self.lap_counter_text.set_text(format_time(current_race_time))
                    
                    # --- HUD UPDATE ---
                    # Update HUD with Selected Driver's Telemetry
                    if hud_driver_name in self.telemetry_data:
                        l_data = self.telemetry_data[hud_driver_name]
                        # Safe Index Check
                        l_len = len(l_data["X"])
                        l_idx = min(idx, l_len - 1)
                        
                        # Speed
                        speed_val = l_data["Speed"][l_idx]
                        self.hud_speed.set_text(f"{int(speed_val)}")
                        
                        # Gear
                        gear_val = int(round(l_data["Gear"][l_idx]))
                        self.hud_gear.set_text(str(gear_val) if gear_val > 0 else "N")
                        self.hud_gear.set_color(l_data["Color"]) # Match driver color
                        
                        # RPM
                        rpm_val = int(l_data["RPM"][l_idx])
                        self.hud_rpm.set_text(f"{rpm_val} RPM")
                        
                        # Throttle (0-100)
                        if self.hud_thr_bar:
                            thr_val = l_data["Throttle"][l_idx]
                            # Width relative to AXES (0.80 max width set in setup_plot)
                            self.hud_thr_bar.set_width(0.80 * (thr_val / 100.0))
                            if self.hud_thr_text:
                                self.hud_thr_text.set_text(f"{int(thr_val)}%")
                        
                        # Brake (0-100 or bool)
                        if self.hud_brk_bar:
                            brk_val = l_data["Brake"][l_idx]
                            if brk_val <= 1.05: brk_val *= 100.0
                            if brk_val > 100: brk_val = 100 
                            self.hud_brk_bar.set_width(0.80 * (brk_val / 100.0))
                            if self.hud_brk_text:
                                self.hud_brk_text.set_text(f"{int(brk_val)}%")
                        
                        # Label Update
                        if self.selected_driver:
                             self.hud_driver.set_text(f"{hud_driver_name}")
                        else:
                             self.hud_driver.set_text(f"{hud_driver_name} (LEADER)")
                        self.hud_driver.set_color(l_data["Color"])

            except Exception as e:
                print(f"Update Error: {e}")
                pass
            
            # DEBUG: Print Lap Number on Change
            if 'last_logged_lap' not in self.__dict__:
                 self.last_logged_lap = -1
            
            if current_lap_num != self.last_logged_lap:
                 print(f"[DEBUG] Current Lap: {current_lap_num} / {total_laps} (Time: {current_race_time:.2f}s)")
                 self.last_logged_lap = current_lap_num
            
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

            # Update Leaderboard (Every frame for smoothness)
            if frame_idx % 1 == 0:
                # SORTING LOGIC: Total Distance
                leaderboard_data.sort(key=lambda x: x[1], reverse=True)
                
                # STORE REFERENCE FOR CLICK HANDLER
                self.leaderboard_data_ref = leaderboard_data
                
                active_leaders = [d for d in leaderboard_data if not d[2]]
                leader_dist = active_leaders[0][1] if active_leaders else 0
                
                # Determine HUD Driver (Selected or Leader)
                # No longer needed to define here as we defined it above, but we keep it for reference if needed
                # hud_driver_name = self.selected_driver if self.selected_driver else leader_driver
                
                for i, (row_frame, pos_lbl, drv_lbl, tyre_lbl, gap_lbl) in enumerate(self.lb_rows):
                    if i < len(leaderboard_data):
                        row_frame.pack(fill="x", pady=2)
                        drv, dist, is_dnf = leaderboard_data[i]
                        
                        # HIGHLIGHT SELECTION
                        if drv == self.selected_driver:
                             # Modern Highlight: Dark Grey with Green Accent/Border
                             row_frame.configure(fg_color="#333333", border_width=1, border_color="#2CC985")
                             pos_lbl.configure(text_color="#2CC985")
                             drv_lbl.configure(text_color="#2CC985", font=("Roboto", 12, "bold"))
                        else:
                             # Default (Dark Mode)
                             row_frame.configure(fg_color="#2b2b2b", border_width=0) 
                             pos_lbl.configure(text_color="white")
                             drv_lbl.configure(text_color="white", font=("Roboto", 12))

                        drv_lbl.configure(text=drv)
                        
                        # Tyre Logic
                        tyre_text = "?"
                        tyre_color = "white"
                        try:
                            # Find current lap for driver
                            d_laps = self.telemetry_data[drv]['Laps']
                            
                            compound = "Unknown"
                            # Filter laps that have started using pre-calculated float time
                            # current_race_time is already normalized
                            
                            if 'NormLapStartTime' in d_laps.columns:
                                # Add 1s buffer to current time to ensure we catch the lap start immediately
                                # matching visual pit exit better and preventing lag
                                started_laps = d_laps[d_laps['NormLapStartTime'] <= (current_race_time + 1.0)]
                            else:
                                # Fallback (shouldn't happen if setup correct)
                                current_session_time_dt = pd.Timedelta(seconds=(current_race_time + self.global_start_time))
                                started_laps = d_laps[d_laps['LapStartTime'] <= current_session_time_dt]
                            
                            if not started_laps.empty:
                                # Get the very latest lap
                                current_lap = started_laps.iloc[-1]
                                compound = str(current_lap['Compound']).strip().upper()
                            else:
                                # Fallback to first lap if race hasn't started for them
                                if not d_laps.empty:
                                    current_lap = d_laps.iloc[0]
                                    compound = str(current_lap['Compound']).strip().upper()

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
                            elif 'HYPERSOFT' in compound: # Legacy support
                                tyre_text = "Hyp"
                                tyre_color = "#FFB6C1" # Pink
                            elif 'ULTRASOFT' in compound: # Legacy support
                                tyre_text = "Utr"
                                tyre_color = "#800080" # Purple
                            elif 'SUPERSOFT' in compound: # Legacy support
                                tyre_text = "Sup"
                                tyre_color = "#FF0000" # Red
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

            
            # --- RETURN UPDATED ARTISTS FOR BLIT ---
            artists = list(self.driver_dots.values()) + list(self.driver_labels.values()) + [self.lap_counter_text]
            
            # Append HUD Artists (if they exist)
            if self.hud_speed: artists.append(self.hud_speed)
            if self.hud_gear: artists.append(self.hud_gear)
            if self.hud_rpm: artists.append(self.hud_rpm)
            if self.hud_driver: artists.append(self.hud_driver)
            if self.hud_unit: artists.append(self.hud_unit)
            if self.hud_thr_bar: artists.append(self.hud_thr_bar)
            if self.hud_brk_bar: artists.append(self.hud_brk_bar)
            if self.hud_thr_text: artists.append(self.hud_thr_text)
            if self.hud_brk_text: artists.append(self.hud_brk_text)
            
            return artists

        self.anim = FuncAnimation(self.fig, update, frames=(total_frames - start_frame) // step_frames, 
                                  interval=base_interval, blit=True, repeat=False)
        self.canvas.draw()

if __name__ == "__main__":
    app = RaceReplayerApp()
    app.mainloop()