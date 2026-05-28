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
import subprocess # Added for Sidecar
import json # Added for Sidecar

# Configure output for standard text
if sys.platform.startswith("win") and sys.stdout is not None:
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

from core.utils import (
    format_time, get_driver_color, cache_event_name_from_folder,
    get_cached_events_for_year, get_cached_years_with_races, build_driver_label
)
from core.data_loader import load_race_data, process_telemetry


class RaceReplayerApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("F1 Race Replayer 2025 - Modern Edition")
        self.geometry("1400x900")
        def maximize_window():
            import sys, os
            if sys.platform.startswith("linux") and os.environ.get("XDG_SESSION_TYPE", "").lower() == "wayland":
                try:
                    w = self.winfo_screenwidth()
                    h = self.winfo_screenheight()
                    self.geometry(f"{w}x{h}+0+0")
                except Exception:
                    pass
            else:
                try:
                    self.state('zoomed')
                except Exception:
                    try:
                        self.attributes('-zoomed', True)
                    except Exception:
                        pass
        self.after(1, maximize_window)
        self.protocol("WM_DELETE_WINDOW", self.on_close)
        
        # State for Selection
        self.selected_driver = None
        self.leaderboard_data_ref = [] # To store current frame's leaderboard order
        self.offline_mode_detected = False
        
        # Animation State
        self.anim = None



        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)

        # -- 1. Header / Controls --
        self.control_frame = ctk.CTkFrame(self, height=60, corner_radius=0)
        self.control_frame.grid(row=0, column=0, columnspan=3, sticky="ew", padx=0, pady=0)
        
        ctk.CTkLabel(self.control_frame, text="Year:", font=("Roboto", 14)).pack(side="left", padx=(20, 5))
        # Dynamic Year Selection
        self.years_list = [str(y) for y in range(2025, 2017, -1)]
        self.year_var = ctk.StringVar(value="2023")
        self.year_combo = ctk.CTkComboBox(self.control_frame, values=self.years_list, variable=self.year_var, width=80, command=self.on_year_changed)
        self.year_combo.pack(side="left", padx=5)
        
        ctk.CTkLabel(self.control_frame, text="Circuit:", font=("Roboto", 14)).pack(side="left", padx=(20, 5))
        
        # Dynamic Circuit Selection
        self.circuit_var = ctk.StringVar(value="Loading...")
        self.circuit_combo = ctk.CTkComboBox(self.control_frame, values=["Loading..."], variable=self.circuit_var, width=180)
        self.circuit_combo.pack(side="left", padx=5)

        ctk.CTkLabel(self.control_frame, text="Driver:", font=("Roboto", 14)).pack(side="left", padx=(20, 5))
        self.driver_var = ctk.StringVar(value="LEADER")
        self.driver_combo = ctk.CTkComboBox(self.control_frame, values=["LEADER"], variable=self.driver_var, width=180, command=self.on_driver_changed)
        self.driver_combo.pack(side="left", padx=5)
        
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
        self.switch_3d.pack(side="right", padx=20)

        # -- 1.1 Session Info Panel (Weather & Status) --
        self.session_info_bar = ctk.CTkFrame(self, height=34, corner_radius=0)
        self.session_info_bar.grid(row=1, column=0, columnspan=3, sticky="ew", padx=0, pady=0)

        self.session_info_frame = ctk.CTkFrame(self.session_info_bar, fg_color="transparent")
        self.session_info_frame.pack(side="right", padx=20)
        
        self.weather_lbl = ctk.CTkLabel(self.session_info_frame, text="Weather: --°C | --%", font=("Mono", 12))
        self.weather_lbl.pack(side="left", padx=(0, 18))
        
        self.track_status_lbl = ctk.CTkLabel(self.session_info_frame, text="TRACK: GREEN", font=("Roboto", 12, "bold"), text_color="#2CC985")
        self.track_status_lbl.pack(side="left")

        # -- 2. Leaderboard (Left Side) --
        self.leaderboard_frame = ctk.CTkFrame(self, width=400, corner_radius=0) # Widened
        self.leaderboard_frame.grid(row=2, column=0, sticky="nsew", padx=0, pady=0)
        
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
        self.map_frame.grid(row=2, column=1, columnspan=2, sticky="nsew")
        
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
        
        # Sidecar State
        self.sidecar_process = None
        self.using_sidecar = False
        self.sidecar_listener_thread = None
        
        # UI State Cache (Prevent Redundant Updates)
        self.row_states = {} 

        
        # Resize Handling
        self.resize_timer = None
        self.was_playing = False
        self.map_frame.bind("<Configure>", self.on_resize)

        # Bind scroll event for 3D zoom
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        
        # --- SMOOTH 3D CAMERA ---
        # State
        self.cam_azim = -60.0
        self.cam_elev = 30.0
        self.target_azim = -60.0
        self.target_elev = 30.0
        self.last_mouse_x = 0
        self.last_mouse_y = 0
        self.is_dragging = False
        
        # Bindings (On the Canvas Widget)
        # Note: We bind to the Tkinter widget to bypass Matplotlib's default handler if needed
        tk_canvas = self.canvas.get_tk_widget()
        tk_canvas.bind("<ButtonPress-1>", self.on_cam_press)
        tk_canvas.bind("<B1-Motion>", self.on_cam_drag)
        tk_canvas.bind("<ButtonRelease-1>", self.on_cam_release)
        
        # Start Physics Loop
        self.update_camera_physics()

        # Initial Circuit Load
        self.on_year_changed(self.year_var.get())

    def on_year_changed(self, new_year):
        """Triggered when user selects a new year from the dropdown."""
        self.circuit_combo.configure(state="disabled", values=["Loading circuits..."])
        self.circuit_var.set("Loading circuits...")
        
        # Launch background thread to fetch circuits
        threading.Thread(target=self.fetch_circuits, args=(int(new_year),), daemon=True).start()

    def fetch_circuits(self, year):
        """Runs in background thread to avoid freezing UI."""
        try:
            # Quick check if cache needs to be explicitly enabled here though it's global
            schedule = fastf1.get_event_schedule(year)
            # Filter out testing sessions robustly using EventFormat
            valid_events = schedule[schedule['EventFormat'] != 'testing']
            circuits = valid_events['EventName'].tolist()
            
            # Update UI on main thread
            self.after(0, self._update_circuit_dropdown, circuits)
        except Exception as e:
            print(f"Error fetching schedule for {year}: {e}")
            self.offline_mode_detected = True
            cached_events = get_cached_events_for_year(year)
            if cached_events:
                cached_years = get_cached_years_with_races()
                if cached_years:
                    self.after(0, self._update_year_dropdown_offline, cached_years, str(year))
                self.after(0, self._update_circuit_dropdown, cached_events)
            else:
                self.after(0, self._update_circuit_dropdown, ["Error loading"])

    def _update_year_dropdown_offline(self, cached_years, requested_year):
        """Switch the year dropdown to cache-backed years when offline."""
        self.year_combo.configure(values=cached_years)

        if requested_year in cached_years:
            self.year_var.set(requested_year)
        elif cached_years:
            self.year_var.set(cached_years[0])

    def _update_circuit_dropdown(self, circuits):
        """Runs on main thread to update Tkinter widgets safely."""
        self.circuit_combo.configure(state="normal", values=circuits)
        if circuits and circuits[0] != "Error loading":
            self.circuit_var.set(circuits[-1]) # Default to last race of the year (usually most interesting)
        else:
            self.circuit_var.set("Not Available")

    def on_driver_changed(self, choice):
        """Update the selected driver from the driver dropdown."""
        if choice == "LEADER" or choice == "Not Available":
            self.selected_driver = None
            self.driver_var.set("LEADER")
            return

        driver_code = choice.split(" - ", 1)[0].strip()
        if driver_code in self.telemetry_data:
            self.selected_driver = driver_code
        else:
            self.selected_driver = None
            self.driver_var.set("LEADER")

    def _update_driver_dropdown(self):
        """Populate driver dropdown after a session has been loaded."""
        driver_codes = list(self.telemetry_data.keys()) if self.telemetry_data else []
        if not driver_codes:
            self.driver_combo.configure(state="disabled", values=["LEADER"])
            self.driver_var.set("LEADER")
            self.selected_driver = None
            return

        driver_codes.sort()
        driver_labels = ["LEADER"] + [build_driver_label(driver_code, self.telemetry_data) for driver_code in driver_codes]
        self.driver_combo.configure(state="normal", values=driver_labels)

        if self.selected_driver in self.telemetry_data:
            selected_label = build_driver_label(self.selected_driver, self.telemetry_data)
            self.driver_var.set(selected_label)
        else:
            self.selected_driver = None
            self.driver_var.set("LEADER")

    def on_cam_press(self, event):
        if not self.is_3d_mode.get(): return
        self.is_dragging = True
        self.last_mouse_x = event.x
        self.last_mouse_y = event.y

    def on_cam_drag(self, event):
        if not self.is_3d_mode.get() or not self.is_dragging: return
        
        dx = event.x - self.last_mouse_x
        dy = event.y - self.last_mouse_y
        
        # Sensitivity
        factor = 0.5
        
        # Update TARGET (not current) to create "pull" effect
        self.target_azim -= dx * factor
        self.target_elev += dy * factor
        
        # Clamp Elev
        self.target_elev = max(0, min(90, self.target_elev))
        
        self.last_mouse_x = event.x
        self.last_mouse_y = event.y

    def on_cam_release(self, event):
        self.is_dragging = False

    def update_camera_physics(self):
        """ Smoothly interpolate camera to target position """
        if self.is_3d_mode.get() and self.ax:
             # Inertia / Damping Factor (Lower = Smoother/Slower)
             lerp_factor = 0.1
             
             # Calculate Diff
             diff_azim = self.target_azim - self.cam_azim
             diff_elev = self.target_elev - self.cam_elev
             
             # Apply if significant movement needed (prevent micro-jitter)
             if abs(diff_azim) > 0.01 or abs(diff_elev) > 0.01:
                 self.cam_azim += diff_azim * lerp_factor
                 self.cam_elev += diff_elev * lerp_factor
                 
                 try:
                     self.ax.view_init(elev=self.cam_elev, azim=self.cam_azim)
                     
                     # Check if we should draw manually
                     # If playing, FuncAnimation is already drawing loop. Drawing here triggers DOUBLE DRAW (Lag).
                     if self.anim and not self.is_paused and not self.is_finished:
                          pass # Do nothing, let animation frame pick up the new angle
                     else:
                          self.canvas.draw_idle() 
                 except: pass

        # Run loop at ~60FPS (16ms)
        self.after(16, self.update_camera_physics)

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

        if self.selected_driver and self.selected_driver in self.telemetry_data:
            self.driver_var.set(build_driver_label(self.selected_driver, self.telemetry_data))
        else:
            self.driver_var.set("LEADER")
            
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
        self.stop_sidecar()
        self.quit()
        self.destroy()
        sys.exit()

    def stop_sidecar(self):
        if self.sidecar_process:
            try:
                # Send gentle exit command
                cmd = json.dumps({"cmd": "exit"}) + "\n"
                self.sidecar_process.stdin.write(cmd)
                self.sidecar_process.stdin.flush()
            except: pass
            
            try:
                self.sidecar_process.terminate()
                self.sidecar_process.wait(timeout=1)
            except: 
                try: self.sidecar_process.kill()
                except: pass
        self.sidecar_process = None
        self.using_sidecar = False
        self.sidecar_listener_thread = None

    def start_replay(self):
        try:
            year = int(self.year_var.get())
            circuit = self.circuit_var.get()
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

        self._update_driver_dropdown()
             
        self.setup_plot()
        self.start_animation(start_frame=int(self.race_start_offset / 0.2))
        self.pause_btn.configure(state="normal", text="PAUSE", fg_color="#E10600")
        self.is_paused = False
        self.is_finished = False
        
        year = self.year_var.get()
        circuit = self.circuit_var.get()
        self.status_lbl.configure(text=f"Replaying: {year} {circuit}", text_color="#2CC985")

                                    
    def toggle_3d_mode(self):
        use_3d = self.is_3d_mode.get()
        if use_3d:
            # Try launching Sidecar
            success = self.start_sidecar()
            if not success:
                 print("Sidecar failed. Falling back to internal Matplotlib 3D.")
                 self.using_sidecar = False
                 # We stay in "3D Mode" but use matplotlib fallback
            else:
                 print("Sidecar launched successfully.")
        else:
            # Stop Sidecar
            self.stop_sidecar()
            
        # Restart plot with new mode if data exists
        if self.telemetry_data:
            self.setup_plot()
            if not self.is_paused and self.anim and self.current_frame < len(self.common_time):
                 self.start_animation(start_frame=self.current_frame)
            else:
                 # Just redraw current frame static
                 pass

    def start_sidecar(self):
        """ Attempts to launch the VisPy Viewer. Returns True if successful. """
        if self.sidecar_process is not None:
             return True # Already running
             
        try:
             # Look for the script in the same directory
             viewer_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "f1_3d_viewer.py")
             # In a frozen PyInstaller bundle, the .py file doesn't physically exist on disk
             if not getattr(sys, 'frozen', False):
                 if not os.path.exists(viewer_script):
                     print(f"Error: Could not find {viewer_script}")
                     return False
                 
             # Launch Subprocess
             # Use the same python executable running this script
             python_exe = sys.executable
             
             kwargs = {
                 'stdin': subprocess.PIPE,
                 'stdout': subprocess.PIPE,
                 'stderr': subprocess.PIPE,
                 'text': True
             }
             if sys.platform == "win32":
                 kwargs['creationflags'] = subprocess.CREATE_NO_WINDOW
             
             self.sidecar_process = subprocess.Popen(
                 [python_exe, viewer_script],
                 **kwargs
             )
             
             # Wait for READY signal or Error
             # Read a few lines to check
             ready = False
             for _ in range(10): 
                 line = self.sidecar_process.stdout.readline()
                 if "CRITICAL" in line or "Error" in line:
                     print(f"Sidecar Init Error: {line}")
                     return False
                 if "VISPY_READY" in line:
                     ready = True
                     break
                     
             if not ready:
                 print("Sidecar did not signal readiness.")
                 self.sidecar_process.terminate()
                 self.sidecar_process = None
                 return False
                 
             self.using_sidecar = True
             
             # --- Send Initial Track Data ---
             if self.telemetry_data:
                 ref_driver = list(self.telemetry_data.keys())[0]
                 ref_x = self.telemetry_data[ref_driver]["X"]
                 ref_y = self.telemetry_data[ref_driver]["Y"]
                 ref_z = self.telemetry_data[ref_driver]["Z"]
                 mask = ~np.isnan(ref_x)
                 valid_x = ref_x[mask]
                 valid_y = ref_y[mask]
                 valid_z = ref_z[mask]
                 
                 step_mesh = 1 # Send FULL resolution track data to VisPy for perfectly smooth lines
                 track_data = {
                     "cmd": "init_track",
                     "x": valid_x[::step_mesh].tolist(),
                     "y": valid_y[::step_mesh].tolist(),
                     "z": valid_z[::step_mesh].tolist()
                 }
                 
                 # Bounds
                 min_x, max_x = np.min(valid_x) - 100, np.max(valid_x) + 100
                 min_y, max_y = np.min(valid_y) - 100, np.max(valid_y) + 100
                 min_z, max_z = np.min(valid_z) - 20, np.max(valid_z) + 20
                 track_data["bounds"] = (min_x, max_x, min_y, max_y, min_z, max_z)
                 
                 try:
                     cmd_str = json.dumps(track_data) + "\n"
                     self.sidecar_process.stdin.write(cmd_str)
                     self.sidecar_process.stdin.flush()
                 except Exception as e:
                     print(f"Failed to send init track: {e}")
                     self.using_sidecar = False
                     return False
             
             # --- Start IPC Listener Thread ---
             self.sidecar_listener_thread = threading.Thread(target=self.listen_to_sidecar, daemon=True)
             self.sidecar_listener_thread.start()
             
             return True
             
        except Exception as e:
             print(f"Exception launching Sidecar: {e}")
             self.sidecar_process = None
             return False

    def listen_to_sidecar(self):
        """ Runs in background thread to catch commands (Pause, Speed) from VisPy """
        if not self.sidecar_process or not self.sidecar_process.stdout: return
        
        try:
            for line in iter(self.sidecar_process.stdout.readline, ''):
                if not line or not self.using_sidecar: break
                line = line.strip()
                if not line: continue
                
                try:
                    data = json.loads(line)
                    cmd = data.get("cmd")
                    
                    if cmd == "input":
                        action = data.get("action")
                        if action == "toggle_pause":
                            self.after(0, self.toggle_pause)
                        elif action == "speed_up":
                            self.after(0, lambda: self.adjust_speed(1))
                        elif action == "speed_down":
                            self.after(0, lambda: self.adjust_speed(-1))
                except json.JSONDecodeError:
                    pass # Ignore non-json print statements from vispy
        except Exception as e:
            print(f"Sidecar Listener Error: {e}")


    def adjust_speed(self, direction):
        """ Helper to cycle speed up or down """
        speeds = ["1x", "5x", "10x", "20x", "50x"]
        current = self.speed_var.get()
        try:
             idx = speeds.index(current)
             new_idx = max(0, min(len(speeds) - 1, idx + direction))
             self.speed_combo.set(speeds[new_idx])
             self.change_speed(speeds[new_idx])
        except ValueError:
             pass

    def setup_plot(self):
        # safely clear
        self.fig.clf() # Clear figure to reset axes (important for 2D vs 3D)
        
        # Check Mode
        use_3d = self.is_3d_mode.get()
        
        # --- 1. SPLIT FIGURES (Map vs HUD) ---
        # Map: Left 75%
        # HUD: Right 25%
        
        if use_3d and not self.using_sidecar:
            # FALLBACK to internal Matplotlib 3D
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
        
        if use_3d and not self.using_sidecar:
             # --- INTERNAL 3D CURTAIN EFFECT (Fallback Only) ---
             try:
                 # 1. Prepare Data & Downsample for Performance
                 # Take every 5th point for the main line, every 10th for the mesh
                 # This reduces vertex count by ~90% for smoother rotation
                 
                 step_line = 5
                 step_mesh = 5 # Medium density 
                 
                 valid_x, valid_y, valid_z = ref_x[mask], ref_y[mask], ref_z[mask]
                 
                 x_mesh_data = valid_x[::step_mesh]
                 y_mesh_data = valid_y[::step_mesh]
                 z_mesh_data = valid_z[::step_mesh]
                 
                 # Base Level
                 z_min = np.min(valid_z) - 20 
                 
                 
                 # 2. Optimized "Curtain" (Vertical Lines via NaN separators)
                 # Structure: [x1, x1, nan, x2, x2, nan, ...]
                 # This draws vertical lines from Track to Base efficiently as one Line object.
                 
                 # Create interleaved arrays
                 # We need to stack (Top, Base, Nan)
                 
                 # Top Points
                 xt = x_mesh_data
                 yt = y_mesh_data
                 zt = z_mesh_data
                 
                 # Base Points
                 xb = x_mesh_data
                 yb = y_mesh_data
                 zb = np.full_like(z_mesh_data, z_min)
                 
                 # NaNs
                 xn = np.full_like(xt, np.nan)
                 yn = np.full_like(yt, np.nan)
                 zn = np.full_like(zt, np.nan)
                 
                 # Stack and Flatten: (3, N) -> (3N,)
                 # We want column-wise stacking: [x1, x2...], [x1, x2...], [nan, nan...]
                 # Flatten 'F' gives: x1, x1, nan, x2, x2, nan...
                 
                 X_lines = np.vstack([xt, xb, xn]).flatten(order='F')
                 Y_lines = np.vstack([yt, yb, yn]).flatten(order='F')
                 Z_lines = np.vstack([zt, zb, zn]).flatten(order='F')
                 
                 # 3. Plot Single Line Collection
                 # CONTRAST FIX: #1C1C1C (Very subtle, just above background #121212)
                 self.ax.plot(X_lines, Y_lines, Z_lines, 
                              color='#1C1C1C', alpha=0.4, linewidth=1)
                 
                 # 4. Plot Downsampled Track Line
                 self.ax.plot(valid_x[::step_line], valid_y[::step_line], valid_z[::step_line], 
                              color='#333333', linewidth=4, alpha=0.5)

             except Exception as e:
                 print(f"3D Curtain Error: {e}")
                 # Fallback
                 self.ax.plot(ref_x[mask], ref_y[mask], ref_z[mask], color='#333333', linewidth=4, alpha=0.5)
        else:
             # 2D Mode OR Sidecar Mode (Draw 2D track on main window)
             self.ax.plot(ref_x[mask], ref_y[mask], color='#333333', linewidth=6, alpha=0.5)
        
        # Lap Counter (Top Left of Map)
        if use_3d and not self.using_sidecar:
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
            if use_3d and not self.using_sidecar:
                 # Internal 3D Plot
                 dot, = self.ax.plot([], [], [], marker='o', color=color, markersize=6, markeredgecolor='black', markeredgewidth=1)
            else:
                 # 2D Plot
                 dot, = self.ax.plot([], [], marker='o', color=color, markersize=8, markeredgecolor='black', markeredgewidth=1)
            
            self.driver_dots[driver] = dot
            
            # Text Setup
            if use_3d and not self.using_sidecar:
                 text = self.ax.text(0, 0, 0, driver, color=color, fontsize=8, fontweight='bold', clip_on=True)
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
        min_x, max_x = np.min(valid_x) - 100, np.max(valid_x) + 100
        min_y, max_y = np.min(valid_y) - 100, np.max(valid_y) + 100
        min_z, max_z = np.min(valid_z) - 20, np.max(valid_z) + 20
        
        self.plot_bounds = (min_x, max_x, min_y, max_y, min_z, max_z)
        self.zoom_level = 1.0
        
        self.plot_bounds = (min_x, max_x, min_y, max_y, min_z, max_z)
        self.zoom_level = 1.0
        
        if use_3d and not self.using_sidecar:
             self.ax.set_zlim(min_z, max_z)
             # Set view
             self.ax.view_init(elev=self.cam_elev, azim=self.cam_azim) # Use stored camera pos
             self.ax.dist = 8 # Move camera closer (Default ~10) to fill screen
             
             # DYNAMIC ASPECT RATIO FIX:
             # Previous (1, 1, 0.2) forced a square box even if track was long/thin.
             # Now we calculate the ratio so it fills the screen.
             x_range = max_x - min_x
             y_range = max_y - min_y
             
             # Maximize relative to the longest side
             max_range = max(x_range, y_range)
             x_ratio = x_range / max_range
             y_ratio = y_range / max_range
             
             print(f"[DEBUG] Aspect Ratio: L={x_range:.0f} W={y_range:.0f} -> ({x_ratio:.2f}, {y_ratio:.2f})")
             
             # Flatten Z (0.2 is generic, maybe scale it too?)
             # 0.2 is usually fine for track elevation vs length
             self.ax.set_box_aspect((x_ratio, y_ratio, 0.2)) 
             
             self.ax.set_axis_off() # Cleaner look
             
             # DISABLE DEFAULT NAV to let our custom physics handler work
             self.ax.finish_resize = lambda: None # Monkeypatch to prevent resets if needed, but set_navigate is key
             # Let's try explicit disable:
             try: self.fig.canvas.toolbar.set_cursor(1) # Reset cursor?
             except: pass
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

            leaderboard_data = [] # moved further up, need to fix
            current_race_time = self.common_time[idx]
            sidecar_frame_data = {}

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
                    sidecar_frame_data[driver] = {"visible": False, "x": 0, "y": 0, "z": 0, "color": d_data["Color"]}
                else:
                    dot.set_visible(True)
                    self.driver_labels[driver].set_visible(True)
                    sidecar_frame_data[driver] = {"visible": True, "x": float(x), "y": float(y), "z": float(z), "color": d_data["Color"]}
                    
                    if self.is_3d_mode.get() and not self.using_sidecar:
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
                    
                    # Update Sidecar Data (Cast is_dnf to pure bool for JSON)
                    if self.using_sidecar:
                         sidecar_frame_data[driver]["is_dnf"] = bool(is_dnf)

            # --- STREAM TO SIDECAR ---
            if self.using_sidecar and self.sidecar_process:
                 # Check if the process is still alive. If user clicked the 'X' on sidecar window,
                 # poll() will return the exit code (not None).
                 if self.sidecar_process.poll() is not None:
                      print("Sidecar window closed by user. Toggling 3D mode OFF.")
                      # Stop sending data
                      self.using_sidecar = False
                      self.sidecar_process = None
                      # We must schedule the UI update back on the main thread
                      self.after(0, lambda: self.switch_3d.deselect()) # Visual update
                      self.after(0, self.toggle_3d_mode) # Internal logic update
                 else:
                      try:
                          # Create UI State Payload
                          minutes = int(current_race_time // 60)
                          seconds = int(current_race_time % 60)
                          ms = int((current_race_time % 1) * 1000)
                          lap_str = f"{minutes:02d}:{seconds:02d}.{ms:03d}"
                          
                          finish_text = ""
                          if self.is_finished: finish_text = "RACE FINISHED"
                          
                          ui_payload = {
                              "speed": self.speed_var.get(),
                              "lap_time": lap_str,
                              "is_paused": self.is_paused,
                              "finish_status": finish_text
                          }
                          
                          msg = json.dumps({"cmd": "update_frame", "drivers": sidecar_frame_data, "ui": ui_payload}) + "\n"
                          self.sidecar_process.stdin.write(msg)
                          self.sidecar_process.stdin.flush()
                      except Exception as e:
                          print(f"Failed to stream frame to sidecar: {e}")
                          # Disable sidecar if pipe breaks
                          self.using_sidecar = False


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
            # Update Session Info (Weather & Status)
            # Optimized: Run every frame but use fast binary search
            try:
                current_session_time = pd.Timedelta(seconds=(current_race_time + self.global_start_time))
                
                # Weather (Keep slower update for weather as it changes slowly)
                if frame_idx % 60 == 0 and self.weather_data is not None:
                    # Find closest weather row before current time
                    # bisect_right returns insertion point to maintain order, -1 gives the last valid entry
                    w_idx = self.weather_data['Time'].searchsorted(current_session_time, side='right') - 1
                    if w_idx >= 0:
                        w_row = self.weather_data.iloc[w_idx]
                        air_temp = w_row['AirTemp']
                        humidity = w_row['Humidity']
                        self.weather_lbl.configure(text=f"Air: {air_temp}°C | Hum: {humidity}%")
                
                # Track Status (Instant Update)
                if self.track_status_data is not None:
                     # Binary Search for fast lookup
                    ts_idx = self.track_status_data['Time'].searchsorted(current_session_time, side='right') - 1
                    
                    if ts_idx >= 0:
                        ts_row = self.track_status_data.iloc[ts_idx]
                        status_code = ts_row['Status']
                        
                        status_text = "GREEN"
                        status_color = "#2CC985"
                        
                        if status_code == '1': 
                            status_text = "GREEN"
                            status_color = "#2CC985"
                        elif status_code == '2':
                            status_text = "YELLOW"
                            status_color = "#FFFF00" # Explicit Yellow
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
            except Exception as e:
                # print(f"Status Error: {e}")
                pass

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
                        drv, dist, is_dnf = leaderboard_data[i]
                        
                        # Calculate Data
                        is_selected = (drv == self.selected_driver)
                        
                        # Tyre Logic
                        tyre_text = "?"
                        tyre_color = "white"
                        try:
                            # Find current lap for driver
                            d_laps = self.telemetry_data[drv]['Laps']
                            
                            compound = "Unknown"
                            if 'NormLapStartTime' in d_laps.columns:
                                started_laps = d_laps[d_laps['NormLapStartTime'] <= (current_race_time + 1.0)]
                            else:
                                current_session_time_dt = pd.Timedelta(seconds=(current_race_time + self.global_start_time))
                                started_laps = d_laps[d_laps['LapStartTime'] <= current_session_time_dt]
                            
                            if not started_laps.empty:
                                current_lap = started_laps.iloc[-1]
                                compound = str(current_lap['Compound']).strip().upper()
                            elif not d_laps.empty:
                                current_lap = d_laps.iloc[0]
                                compound = str(current_lap['Compound']).strip().upper()

                            if 'SOFT' in compound:
                                tyre_text = "Soft"; tyre_color = "#FF3333"
                            elif 'MEDIUM' in compound:
                                tyre_text = "Medium"; tyre_color = "#FFE000"
                            elif 'HARD' in compound:
                                tyre_text = "Hard"; tyre_color = "white"
                            elif 'INTER' in compound:
                                tyre_text = "Inter"; tyre_color = "#39B54A"
                            elif 'WET' in compound:
                                tyre_text = "Wet"; tyre_color = "#00AEEF"
                            elif 'HYPERSOFT' in compound:
                                tyre_text = "Hyp"; tyre_color = "#FFB6C1"
                            elif 'ULTRASOFT' in compound:
                                tyre_text = "Utr"; tyre_color = "#800080"
                            elif 'SUPERSOFT' in compound:
                                tyre_text = "Sup"; tyre_color = "#FF0000"
                            else:
                                tyre_text = compound[:3] if compound and compound != "NAN" else "?"
                                tyre_color = "gray"
                        except: 
                            tyre_text = "Err"

                        # Gap Logic
                        if is_dnf:
                             gap_text = "DNF"
                             gap_color = "red"
                             main_text_color = "gray"
                        else:
                            main_text_color = "white"
                            if i == 0:
                                gap_text = "LEADER"
                                gap_color = "#2CC985"
                                main_text_color = "#2CC985"
                            else:
                                gap = leader_dist - dist
                                gap_text = f"-{gap:.0f} m"
                                gap_color = "white"
                        
                        if is_selected:
                             main_text_color = "#2CC985"
                        
                        # --- STATE CHECK ---
                        # Create a signature of the visual state
                        new_state = (drv, is_selected, tyre_text, tyre_color, gap_text, gap_color, main_text_color)
                        
                        if self.row_states.get(i) == new_state:
                            continue # SKIP UPDATE IF IDENTICAL
                            
                        # Apply Updates
                        self.row_states[i] = new_state
                        row_frame.pack(fill="x", pady=2) # Ensure packed
                        
                        # 1. Style Container
                        if is_selected:
                             row_frame.configure(fg_color="#333333", border_width=1, border_color="#2CC985")
                             drv_lbl.configure(font=("Roboto", 12, "bold"))
                        else:
                             # Transparent makes it look cleaner and prevents "box" artifacts
                             row_frame.configure(fg_color="transparent", border_width=0) 
                             drv_lbl.configure(font=("Roboto", 12))

                        # 2. Update Text Content
                        drv_lbl.configure(text=drv, text_color=main_text_color)
                        pos_lbl.configure(text_color=main_text_color)
                        
                        tyre_lbl.configure(text=tyre_text, text_color=tyre_color)
                        gap_lbl.configure(text=gap_text, text_color=gap_color)
                        
                    else:
                        row_frame.pack_forget()
                        self.row_states[i] = None # Clear state

            
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

        # Determine Blit Setting dynamically
        raw_3d = self.is_3d_mode.get()
        is_3d = (str(raw_3d) == "1" or str(raw_3d).lower() == "on" or raw_3d is True)
        use_blit = not is_3d
        print(f"[DEBUG] Animation Start: 3D={is_3d}, Blit={use_blit}")

        self.anim = FuncAnimation(self.fig, update, frames=(total_frames - start_frame) // step_frames, 
                                  interval=base_interval, blit=use_blit, repeat=False)
        self.canvas.draw()

if __name__ == "__main__":
    app = RaceReplayerApp()
    app.mainloop()