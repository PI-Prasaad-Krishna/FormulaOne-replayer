import os
import fastf1
import fastf1.plotting

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

def cache_event_name_from_folder(folder_name):
    """Convert a cached FastF1 event folder into a readable event name."""
    if "_" not in folder_name:
        return folder_name.replace("_", " ")

    parts = folder_name.split("_", 1)
    if len(parts) != 2:
        return folder_name.replace("_", " ")

    return parts[1].replace("_", " ")

def get_cached_events_for_year(year):
    """Return cached event names for a given year when online schedules are unavailable."""
    year_dir = os.path.join("f1_cache", str(year))
    if not os.path.isdir(year_dir):
        return []

    events = []
    seen = set()

    try:
        for entry in sorted(os.listdir(year_dir)):
            event_dir = os.path.join(year_dir, entry)
            if not os.path.isdir(event_dir):
                continue

            for session_dir in os.listdir(event_dir):
                if not session_dir.endswith("_Race"):
                    continue

                event_name = cache_event_name_from_folder(entry)
                if event_name not in seen:
                    seen.add(event_name)
                    events.append(event_name)
                break
    except Exception:
        return []

    return events

def get_cached_years_with_races():
    """Return cached years that contain at least one race session."""
    cache_root = "f1_cache"
    if not os.path.isdir(cache_root):
        return []

    years = []
    try:
        for entry in sorted(os.listdir(cache_root), reverse=True):
            if not entry.isdigit():
                continue

            year_dir = os.path.join(cache_root, entry)
            if not os.path.isdir(year_dir):
                continue

            has_race = False
            for event_entry in os.listdir(year_dir):
                event_dir = os.path.join(year_dir, event_entry)
                if not os.path.isdir(event_dir):
                    continue

                for session_dir in os.listdir(event_dir):
                    if session_dir.endswith("_Race"):
                        has_race = True
                        break

                if has_race:
                    break

            if has_race:
                years.append(entry)
    except Exception:
        return []

    return years

def build_driver_label(driver_code, telemetry_data):
    data = telemetry_data.get(driver_code, {})
    team = data.get("Team", "")
    if team and team != "Unknown":
        return f"{driver_code} - {team}"
    return driver_code
