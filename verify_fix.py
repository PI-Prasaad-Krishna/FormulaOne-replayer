import fastf1
import pandas as pd
import numpy as np
import os

# Setup Cache
if not os.path.exists("f1_cache"):
    os.makedirs("f1_cache")
fastf1.Cache.enable_cache("f1_cache")

def verify_fix(year, circuit, session_type="R"):
    print(f"Loading {year} {circuit} [{session_type}] session data...")
    try:
        session = fastf1.get_session(year, circuit, session_type)
        session.load(telemetry=True, laps=True, weather=False)
    except Exception as e:
        print(f"Error loading session: {e}")
        return

    drivers = pd.unique(session.laps['Driver'])
    
    # --- 1. Global Start Time Logic ---
    global_start_time = None
    ref_driver = None
    max_laps = 0

    for d in drivers:
        try:
            n = session.laps.pick_driver(d)['LapNumber'].max()
            if n > max_laps:
                max_laps = n
                ref_driver = d
        except: pass
    
    if ref_driver:
        try:
             ref_laps = session.laps.pick_driver(ref_driver)
             lap1 = ref_laps[ref_laps['LapNumber'] == 1]
             if not lap1.empty:
                 global_start_time = lap1.iloc[0]['LapStartTime'].total_seconds()
                 print(f"Global Start Time (Lap 1 Start): {global_start_time}s")
        except Exception as e:
             print(f"Error getting global start: {e}")

    if global_start_time is None:
        print("FAILED to determine global start time")
        return

    # --- 2. Max Race Time Logic ---
    max_race_time = 0
    
    for driver in drivers:
        try:
            laps = session.laps.pick_driver(driver)
            
            # Simulated Telemetry Check
            if not laps.empty:
                tel = laps.get_telemetry()
                if not tel.empty:
                    t_seconds = tel['Time'].dt.total_seconds().to_numpy() - global_start_time
                    if len(t_seconds) > 0:
                         max_race_time = max(max_race_time, t_seconds[-1])
            
            # Robust Lap Check (The Fix)
            if not laps.empty:
                last_lap = laps.iloc[-1]
                if pd.notna(last_lap['LapStartTime']) and pd.notna(last_lap['LapTime']):
                    end_t = (last_lap['LapStartTime'] + last_lap['LapTime']).total_seconds() - global_start_time
                    if end_t > max_race_time:
                        max_race_time = end_t
                        
        except Exception:
            pass
            
    print(f"Calculated Max Race Time: {max_race_time:.2f}s")
    
    # --- 3. Lap Counter Logic ---
    ref_laps = session.laps.pick_driver(ref_driver)
    lap_start_times = ref_laps['LapStartTime'].dt.total_seconds().to_numpy() - global_start_time
    lap_numbers = ref_laps['LapNumber'].to_numpy()
    
    lap_start_times = np.nan_to_num(lap_start_times, nan=0.0)
    sort_idx = np.argsort(lap_start_times)
    lap_start_times = lap_start_times[sort_idx]
    lap_numbers = lap_numbers[sort_idx]
    
    print("\nLap Counter alignment:")
    print(f"Lap 1 Start Time: {lap_start_times[0]:.2f}s (Should be 0.0)")
    print(f"Final Lap ({lap_numbers[-1]}) Start Time: {lap_start_times[-1]:.2f}s")
    print(f"Does Max Race Time cover Final Lap Start? {max_race_time > lap_start_times[-1]}")
    
    # Check what frame count would be
    step_size = 0.2
    total_frames = int(max_race_time / step_size)
    print(f"Total Frames: {total_frames}")

if __name__ == "__main__":
    verify_fix(2023, "Abu Dhabi")
    print("\n---")
    verify_fix(2024, "Monza") # Approximation for 2025
