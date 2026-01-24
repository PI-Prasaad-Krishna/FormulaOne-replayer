import fastf1
import pandas as pd
import numpy as np
import datetime

# Setup cache
fastf1.Cache.enable_cache("f1_cache")

def debug_lap_data(year, circuit):
    print(f"Loading {year} {circuit}...")
    try:
        session = fastf1.get_session(year, circuit, 'R')
        session.load(telemetry=True, laps=True)
    except Exception as e:
        print(f"Load Error: {e}")
        raise e

    with open("debug_output.txt", "w") as f:
        f.write(f"Loading {year} {circuit}...\n")
        
        # 1. Calculate Global Start Time (Min Telemetry Time)
        min_session_time = float('inf')
        drivers = pd.unique(session.laps['Driver'])
        
        f.write("\n--- Telemetry Start Times ---\n")
        for driver in drivers[:5]: # Check first few drivers
            try:
                laps = session.laps.pick_drivers(driver)
                tel = laps.get_telemetry()
                if not tel.empty:
                    t_min = tel['Time'].dt.total_seconds().min()
                    f.write(f"Driver {driver} Min Time: {t_min}\n")
                    if t_min < min_session_time:
                        min_session_time = t_min
            except Exception as e:
                f.write(f"Error {driver}: {e}\n")
                
        f.write(f"\nGlobal Start Time (Min Session Time): {min_session_time}\n")
        
        # 2. Inspect Leader's Laps
        f.write("\n--- Leader Lap Data ---\n")
        # Find leader (most laps)
        max_laps = 0
        leader = drivers[0]
        for d in drivers:
            try:
                n = session.laps.pick_drivers(d)['LapNumber'].max()
                if n > max_laps:
                    max_laps = n
                    leader = d
            except: pass
                
        f.write(f"Likely Leader: {leader} ({max_laps} laps)\n")
        
        laps = session.laps.pick_drivers(leader)
        
        f.write(f"\n{'Lap':<5} | {'LapStartTime':<25} | {'NormStartTime (s)':<20} | {'LapTime':<20}\n")
        f.write("-" * 80 + "\n")
        
        for _, row in laps.iterrows():
            lap_n = row['LapNumber']
            lst = row['LapStartTime']
            lt = row['LapTime']
            
            norm_start = np.nan
            if pd.notna(lst):
                norm_start = lst.total_seconds() - min_session_time
                
            f.write(f"{lap_n:<5} | {str(lst):<25} | {norm_start:<20.2f} | {str(lt):<20}\n")

if __name__ == "__main__":
    try:
        debug_lap_data(2025, "Japan")
    except Exception as e:
        print(f"Failed 2025 Japan: {e}")
        try:
            debug_lap_data(2023, "Abu Dhabi")
        except Exception as e2:
            print(f"Failed 2023 Abu Dhabi: {e2}")
