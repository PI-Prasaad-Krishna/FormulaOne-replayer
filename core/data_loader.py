import fastf1
import pandas as pd
import numpy as np
from .utils import get_driver_color

def load_race_data(year, circuit, session_type="R"):
    print(f"Loading {year} {circuit} [{session_type}] session data...")
    
    def internal_load():
        session = fastf1.get_session(year, circuit, session_type)
        session.load(telemetry=True, laps=True, weather=True) # EXTENDED: Load Weather
        return session
        
    try:
        try:
            return internal_load()
        except Exception as e1:
            print(f"Online load failed ({e1}). Retrying with Offline Mode...")
            
            try:
                fastf1.Cache.offline_mode(enabled=True)
                return internal_load()
            except Exception as e2:
                print(f"Offline load also failed: {e2}")
                raise e2
            finally:
                fastf1.Cache.offline_mode(enabled=False)
                
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
