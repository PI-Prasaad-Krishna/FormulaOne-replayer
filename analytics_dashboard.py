# F1 Analytics Dashboard v1.0 (FastF1 + CustomTkinter)
# -----------------------------------------------------------------------------
# Features:
# - Telemetry Battle (Head-to-Head)
# - Race Strategy Analysis (Tyre History, Lap Times)
# - Track Visualization (Gear Shifts, Speed Map)
# - Driver DNA (Radar Chart)
# - Position Change Chart
# -----------------------------------------------------------------------------

import customtkinter as ctk
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import fastf1
import fastf1.plotting
import pandas as pd
import numpy as np
import threading
import os
from scipy.interpolate import splprep, splev

# Configure appearance
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

class AnalyticsDashboardApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("F1 Analytics Dashboard 2025")
        self.geometry("1400x900")
        
        # Grid Layout
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)
        
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.is_running = True

        # --- 3D Camera State ---
        self.cam_azim = -60.0
        self.cam_elev = 30.0
        self.target_azim = -60.0
        self.target_elev = 30.0
        self.zoom_level = 1.0
        self.last_mouse_x = 0
        self.last_mouse_y = 0
        self.is_dragging = False
        
        # Start Physics Loop
        self.update_camera_physics()

        # --- Sidebar ---
        self.sidebar = ctk.CTkFrame(self, width=250, corner_radius=0)
        self.sidebar.grid(row=0, column=0, sticky="nsew")
        
        ctk.CTkLabel(self.sidebar, text="ANALYTICS", font=("Roboto", 20, "bold"), text_color="#E10600").pack(pady=(30, 20))
        
        # Session Selection
        ctk.CTkLabel(self.sidebar, text="Year:").pack(padx=20, anchor="w")
        self.year_var = ctk.StringVar(value="2023")
        ctk.CTkEntry(self.sidebar, textvariable=self.year_var).pack(padx=20, pady=(0, 10), fill="x")
        
        ctk.CTkLabel(self.sidebar, text="Circuit:").pack(padx=20, anchor="w")
        self.event_var = ctk.StringVar(value="Abu Dhabi")
        ctk.CTkEntry(self.sidebar, textvariable=self.event_var).pack(padx=20, pady=(0, 20), fill="x")
        
        ctk.CTkLabel(self.sidebar, text="PAGES", font=("Roboto", 12, "bold"), text_color="gray").pack(padx=20, anchor="w", pady=(10, 5))
        
        self.btn_telemetry = ctk.CTkButton(self.sidebar, text="Telemetry Battle", command=lambda: self.show_page("telemetry"), fg_color="transparent", border_width=1)
        self.btn_telemetry.pack(pady=10, padx=20, fill="x")
        
        self.btn_strategy = ctk.CTkButton(self.sidebar, text="Race Strategy", command=lambda: self.show_page("strategy"), fg_color="transparent", border_width=1)
        self.btn_strategy.pack(pady=10, padx=20, fill="x")
        
        self.btn_track = ctk.CTkButton(self.sidebar, text="Track Dominance", command=lambda: self.show_page("track"), fg_color="transparent", border_width=1)
        self.btn_track.pack(pady=10, padx=20, fill="x")
        
        self.btn_dna = ctk.CTkButton(self.sidebar, text="Driver DNA", command=lambda: self.show_page("dna"), fg_color="transparent", border_width=1)
        self.btn_dna.pack(pady=10, padx=20, fill="x")

        self.btn_pos = ctk.CTkButton(self.sidebar, text="Position Chart", command=lambda: self.show_page("position"), fg_color="transparent", border_width=1)
        self.btn_pos.pack(pady=10, padx=20, fill="x")

        self.btn_lines = ctk.CTkButton(self.sidebar, text="Race Lines", command=lambda: self.show_page("race_lines"), fg_color="transparent", border_width=1)
        self.btn_lines.pack(pady=10, padx=20, fill="x")

        self.btn_tyre = ctk.CTkButton(self.sidebar, text="Tyre Degradation", command=lambda: self.show_page("tyre"), fg_color="transparent", border_width=1)
        self.btn_tyre.pack(pady=10, padx=20, fill="x")
        
        self.btn_track_info = ctk.CTkButton(self.sidebar, text="Track Info (3D)", command=lambda: self.show_page("track_info"), fg_color="transparent", border_width=1)
        self.btn_track_info.pack(pady=10, padx=20, fill="x")

        # Navigation Buttons Map
        self.nav_buttons = {
            "telemetry": self.btn_telemetry,
            "strategy": self.btn_strategy,
            "track": self.btn_track,
            "dna": self.btn_dna,
            "position": self.btn_pos,
            "race_lines": self.btn_lines,
            "tyre": self.btn_tyre,
            "track_info": self.btn_track_info
        }

        # --- Main Content Area ---
        self.content_area = ctk.CTkFrame(self, corner_radius=0, fg_color="#121212")
        self.content_area.grid(row=0, column=1, sticky="nsew")
        
        # Default Page
        self.show_page("telemetry")

    def show_page(self, page_name):
        # Update styling for sidebar buttons
        for name, btn in self.nav_buttons.items():
            if name == page_name:
                btn.configure(fg_color="#E10600")
            else:
                btn.configure(fg_color="transparent")

        # Clear content area
        for widget in self.content_area.winfo_children():
            widget.destroy()
            
        if page_name == "telemetry":
            self.create_telemetry_page()
        elif page_name == "strategy":
            self.create_strategy_page()
        elif page_name == "track":
            self.create_track_page()
        elif page_name == "dna":
            self.create_dna_page()
        elif page_name == "position":
            self.create_position_page()
        elif page_name == "race_lines":
            self.create_race_lines_page()
        elif page_name == "tyre":
            self.create_tyre_page()
        elif page_name == "track_info":
            self.create_track_info_page()

    def load_session(self, session_type='R'):
        # Helper to load session if not loaded (hardcoded 2023 Abu Dhabi for demo or add inputs)
        # For simplicity, we'll use a dialog or just 2023 Abu Dhabi
        try:
             year = int(self.year_var.get())
             circuit = self.event_var.get()
        except:
             messagebox.showerror("Error", "Invalid Year")
             return None
             
        status_lbl = ctk.CTkLabel(self.sidebar, text=f"Loading {session_type} Data...", text_color="orange")
        status_lbl.pack(side="bottom", pady=10)
        self.update()
        
        try:
            # Enable Cache
            if not os.path.exists("f1_cache"): os.makedirs("f1_cache")
            fastf1.Cache.enable_cache("f1_cache")
            
            # Map full name to identifier if needed, but 'Race'/'Qualifying' -> 'R'/'Q'
            identifier = 'R'
            if session_type.lower().startswith('q'): identifier = 'Q'
            elif session_type.lower().startswith('r'): identifier = 'R'
            else: identifier = session_type # Trace or others
            
            session = fastf1.get_session(year, circuit, identifier)
            session.load(telemetry=True, laps=True, weather=False)
            status_lbl.destroy()
            return session
        except Exception as e:
             status_lbl.configure(text="Error", text_color="red")
             print(f"Load Error: {e}")
             return None

    def get_compound_color_safe(self, compound, session=None):
        """
        Safely retrieve compound color supporting both old (COMPOUND_COLORS) 
        and new (get_compound_color) FastF1 versions.
        """
        # Try new API first
        if hasattr(fastf1.plotting, 'get_compound_color'):
            try:
                return fastf1.plotting.get_compound_color(compound, session=session)
            except:
                pass
        
        # Try old API
        if hasattr(fastf1.plotting, 'COMPOUND_COLORS'):
             return fastf1.plotting.COMPOUND_COLORS.get(compound, "white")
             
        # Fallback manual map for future proofing if both fail
        COMPOUND_MAP = {
            "SOFT": "#da291c",
            "MEDIUM": "#ffd12e",
            "HARD": "#f0f0ec",
            "INTERMEDIATE": "#43d662",
            "WET": "#0096ff",
            "HYPERSOFT": "#ffb3c3",
            "ULTRASOFT": "#b148c9",
            "SUPERSOFT": "#da291c",
            "SOFT": "#ffd12e", # 2018? No, modern
             # Modern Pirelli colors usually:
             # Soft: Red, Medium: Yellow, Hard: WhiteF
            "C1": "#f0f0ec", "C2": "#f0f0ec",
            "C3": "#ffd12e", "C4": "#da291c", "C5": "#da291c"
        }
        # Actually, let's just stick to the basic names as FastF1 usually normalizes them
        return COMPOUND_MAP.get(compound, "white")

    def create_telemetry_page(self):
        ctk.CTkLabel(self.content_area, text="Interactive Telemetry Battle", font=("Roboto", 24, "bold")).pack(pady=10)
        
        controls = ctk.CTkFrame(self.content_area)
        controls.pack(pady=10)
        
        self.d1_var = ctk.StringVar(value="VER")
        self.d2_var = ctk.StringVar(value="HAM")
        self.tel_session_var = ctk.StringVar(value="Race")
        
        ctk.CTkOptionMenu(controls, variable=self.tel_session_var, values=["Race", "Qualifying"], width=100, fg_color="#333", button_color="#444").pack(side="left", padx=10)
        
        ctk.CTkEntry(controls, textvariable=self.d1_var, width=60).pack(side="left", padx=10)
        ctk.CTkLabel(controls, text="VS").pack(side="left", padx=10)
        ctk.CTkEntry(controls, textvariable=self.d2_var, width=60).pack(side="left", padx=10)
        
        ctk.CTkButton(controls, text="Analyze Fastest Lap", command=self.plot_telemetry, fg_color="#E10600").pack(side="left", padx=20)
        self.btn_play = ctk.CTkButton(controls, text="Play Replay", command=self.toggle_telemetry_animation, fg_color="#333", width=100, state="disabled")
        self.btn_play.pack(side="left", padx=10)
        
        self.telemetry_frame = ctk.CTkFrame(self.content_area, fg_color="transparent")
        self.telemetry_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
    def plot_telemetry(self):
        session_type = self.tel_session_var.get()
        session = self.load_session(session_type)
        if not session: return
        
        d1 = self.d1_var.get().upper()
        d2 = self.d2_var.get().upper()
        
        try:
            laps_d1 = session.laps.pick_drivers(d1).pick_fastest()
            laps_d2 = session.laps.pick_drivers(d2).pick_fastest()
            
            if laps_d1 is None or laps_d2 is None:
                messagebox.showerror("Error", "Data not available for one or both drivers.")
                return

            # Unified Data Source (Video/Position/Car Data merged)
            self.tel_d1 = laps_d1.get_telemetry().add_distance()
            self.tel_d2 = laps_d2.get_telemetry().add_distance()
            
            # Position Data is now same as Telemetry Data
            self.tel_d1_pos = self.tel_d1
            self.tel_d2_pos = self.tel_d2
            
            c1 = fastf1.plotting.get_driver_color(d1, session=session)
            c2 = fastf1.plotting.get_driver_color(d2, session=session)
            
            # Enable play button
            self.btn_play.configure(state="normal", fg_color="#333")
            
        except Exception as e:
            messagebox.showerror("Error", f"Analysis error: {e}")
            return
        
        # Plot
        for widget in self.telemetry_frame.winfo_children(): widget.destroy()
        
        # Grid Layout for Plots: Left (Traces), Right (Track)
        fig = plt.figure(figsize=(12, 8), facecolor='#121212')
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[2, 1])
        
        # Speed Trace (Top Left)
        self.ax_speed = fig.add_subplot(gs[0, 0])
        self.ax_speed.set_facecolor('#121212')
        self.ax_speed.plot(self.tel_d1['Distance'], self.tel_d1['Speed'], color=c1, label=d1)
        self.ax_speed.plot(self.tel_d2['Distance'], self.tel_d2['Speed'], color=c2, label=d2)
        self.ax_speed.set_ylabel("Speed (km/h)", color="white")
        self.ax_speed.tick_params(colors="white", labelbottom=False)
        self.ax_speed.grid(True, linestyle='--', alpha=0.2)
        self.ax_speed.legend(facecolor='#333', labelcolor='white')
        self.ax_speed.set_title(f"Telemetry: {d1} vs {d2}", color="white")
        
        # HUD Text (Top Right of Speed Plot)
        self.hud_text_d1 = self.ax_speed.text(0.98, 0.9, f"{d1}: 0 km/h", transform=self.ax_speed.transAxes, 
                                            color=c1, ha='right', fontweight='bold', fontsize=10)
        self.hud_text_d2 = self.ax_speed.text(0.98, 0.8, f"{d2}: 0 km/h", transform=self.ax_speed.transAxes, 
                                            color=c2, ha='right', fontweight='bold', fontsize=10)
        
        # Throttle Trace (Bottom Left)
        self.ax_throt = fig.add_subplot(gs[1, 0], sharex=self.ax_speed)
        self.ax_throt.set_facecolor('#121212')
        self.ax_throt.plot(self.tel_d1['Distance'], self.tel_d1['Throttle'], color=c1, linestyle="--")
        self.ax_throt.plot(self.tel_d2['Distance'], self.tel_d2['Throttle'], color=c2, linestyle="--")
        self.ax_throt.set_ylabel("Throttle %", color="white")
        self.ax_throt.set_xlabel("Distance (m)", color="white")
        self.ax_throt.tick_params(colors="white")
        self.ax_throt.grid(True, linestyle='--', alpha=0.2)
        
        # Track Map (Right Column)
        self.ax_track = fig.add_subplot(gs[:, 1])
        self.ax_track.set_facecolor('#121212')
        self.ax_track.plot(self.tel_d1_pos['X'], self.tel_d1_pos['Y'], color=c1, alpha=0.5, linewidth=1)
        self.ax_track.plot(self.tel_d2_pos['X'], self.tel_d2_pos['Y'], color=c2, alpha=0.5, linewidth=1, linestyle='--')
        self.ax_track.axis('equal')
        self.ax_track.axis('off')
        self.ax_track.set_title("Live Position", color="white")
        
        # Interactive Markers
        self.cursor_speed = self.ax_speed.axvline(x=0, color='white', linestyle=':', alpha=0.8)
        self.cursor_throt = self.ax_throt.axvline(x=0, color='white', linestyle=':', alpha=0.8)
        
        self.marker_d1, = self.ax_track.plot([], [], 'o', color='white', markersize=8, markeredgecolor=c1, markeredgewidth=2)
        self.marker_d2, = self.ax_track.plot([], [], 'o', color='white', markersize=8, markeredgecolor=c2, markeredgewidth=2)
        
        self.tel_canvas = FigureCanvasTkAgg(fig, master=self.telemetry_frame)
        self.tel_canvas.get_tk_widget().pack(fill="both", expand=True)
        self.tel_canvas.draw()
        
        # Connect Hover Event
        # Connect Hover Event
        self.tel_canvas.mpl_connect('motion_notify_event', self.on_hover_telemetry)

    def toggle_telemetry_animation(self):
        if hasattr(self, 'anim_running') and self.anim_running:
            self.anim_running = False
            self.btn_play.configure(text="Play Replay")
        else:
            self.anim_running = True
            self.btn_play.configure(text="Pause Replay")
            self.animate_telemetry()

    def animate_telemetry(self):
        if not hasattr(self, 'anim_running') or not self.anim_running: return
        if not getattr(self, 'is_running', True): return # Stop if window closed
        
        # Determine max distance
        max_dist = self.tel_d1['Distance'].max()
        
        # Current distance state
        if not hasattr(self, 'current_dist'):
            self.current_dist = 0
            
        step = 20 # Meters per frame (adjust for speed)
        
        self.current_dist += step
        if self.current_dist > max_dist:
            self.current_dist = 0 # Loop or stop? Let's loop
            
        self.update_telemetry_cursor(self.current_dist)
        
        # Schedule next frame (50ms = 20fps)
        self.after(50, self.animate_telemetry)

    def update_telemetry_cursor(self, dist):
        if not hasattr(self, 'cursor_speed'): return
        
        # Update Cursors
        self.cursor_speed.set_xdata([dist])
        self.cursor_throt.set_xdata([dist])
        
        try:
            # D1
            idx1 = (np.abs(self.tel_d1_pos['Distance'] - dist)).argmin()
            x1, y1 = self.tel_d1_pos['X'].iloc[idx1], self.tel_d1_pos['Y'].iloc[idx1]
            self.marker_d1.set_data([x1], [y1])
            
            # Update HUD D1
            s1 = self.tel_d1['Speed'].iloc[idx1]
            t1 = self.tel_d1['Throttle'].iloc[idx1]
            self.hud_text_d1.set_text(f"{self.d1_var.get().upper()}: {s1:.0f} km/h | {t1:.0f}%")
            
            # D2
            idx2 = (np.abs(self.tel_d2_pos['Distance'] - dist)).argmin()
            x2, y2 = self.tel_d2_pos['X'].iloc[idx2], self.tel_d2_pos['Y'].iloc[idx2]
            self.marker_d2.set_data([x2], [y2])
            
            # Update HUD D2
            s2 = self.tel_d2['Speed'].iloc[idx2]
            t2 = self.tel_d2['Throttle'].iloc[idx2]
            self.hud_text_d2.set_text(f"{self.d2_var.get().upper()}: {s2:.0f} km/h | {t2:.0f}%")
            
            self.tel_canvas.draw_idle()
        except: pass

    def on_hover_telemetry(self, event):
        if event.inaxes in [self.ax_speed, self.ax_throt]:
            dist = event.xdata
            if dist is None: return
            self.anim_running = False # Stop animation if user interacts
            self.btn_play.configure(text="Play Replay")
            self.update_telemetry_cursor(dist)
            
        elif event.inaxes == self.ax_track:
            mx, my = event.xdata, event.ydata
            if mx is None: return
            self.anim_running = False
            self.btn_play.configure(text="Play Replay")
            
            try:
                dists = np.sqrt((self.tel_d1_pos['X'] - mx)**2 + (self.tel_d1_pos['Y'] - my)**2)
                idx = dists.argmin()
                found_dist = self.tel_d1_pos['Distance'].iloc[idx]
                self.update_telemetry_cursor(found_dist)
            except: pass

    def create_strategy_page(self):
        ctk.CTkLabel(self.content_area, text="Race Strategy Analysis", font=("Roboto", 24, "bold")).pack(pady=10)
        ctk.CTkButton(self.content_area, text="Load Strategy Data", command=self.plot_strategy, fg_color="#E10600").pack(pady=5)
        
        self.strategy_frame = ctk.CTkFrame(self.content_area, fg_color="transparent")
        self.strategy_frame.pack(fill="both", expand=True, padx=10, pady=10)

    def plot_strategy(self):
        session = self.load_session()
        if not session: return
        
        for widget in self.strategy_frame.winfo_children(): widget.destroy()
        
        fig, ax = plt.subplots(2, 1, figsize=(10, 8), facecolor='#121212')
        fig.subplots_adjust(hspace=0.4) # Better spacing
        
        # 1. Lap Time Evolution
        ax[0].set_facecolor('#121212')
        for drv in session.laps['Driver'].unique():
             drv_laps = session.laps.pick_drivers(drv)
             # Filter outliers
             drv_laps = drv_laps.pick_quicklaps(1.07)
             if not drv_laps.empty:
                 color = fastf1.plotting.get_driver_color(drv, session=session)
                 ax[0].scatter(drv_laps['LapNumber'], drv_laps['LapTime'].dt.total_seconds(), 
                               color=color, s=10, alpha=0.6)
                               
        ax[0].set_ylabel("Lap Time (s)", color="white")
        ax[0].set_xlabel("Lap Number", color="white")
        ax[0].tick_params(colors="white")
        ax[0].set_title("Lap Time Evolution", color="white")
        
        # 2. Tyre Strategy (Horizontal Bar)
        ax[1].set_facecolor('#121212')
        
        # Plot stints
        stints = session.laps.groupby(["Driver", "Stint", "Compound"])
        stints = stints.count().reset_index()
        stints = stints[["Driver", "Stint", "Compound", "LapNumber"]]
        stints = stints.rename(columns={"LapNumber": "StintLength"})
        stints = stints.sort_values(by=["Stint"])
        
        # Need start lap of each stint
        # This is complex to robustly calc quickly, but FastF1 has visuals.
        # We will iterate drivers
        drivers = session.drivers
        y_positions = range(len(drivers))
        
        for i, driver in enumerate(drivers):
            d_laps = session.laps.pick_drivers(driver)
            if d_laps.empty: continue
            
            # Identify stints
            # Change in compound or logic
            stints = d_laps[['LapNumber', 'Compound']].groupby((d_laps['Compound'] != d_laps['Compound'].shift()).cumsum())
            
            for _, stint in stints:
                start_lap = stint['LapNumber'].min()
                end_lap = stint['LapNumber'].max()
                length = end_lap - start_lap + 1
                compound = stint['Compound'].iloc[0]
                
                color = self.get_compound_color_safe(compound, session=session)
                
                ax[1].barh(i, length, left=start_lap, color=color, edgecolor='black', height=0.8)
                
        ax[1].set_yticks(y_positions)
        ax[1].set_yticklabels(drivers, color="white", fontsize=8)
        ax[1].set_xlabel("Lap Number", color="white")
        ax[1].tick_params(colors="white")
        ax[1].invert_yaxis()
        ax[1].set_title("Tyre Strategies", color="white")
        
        canvas = FigureCanvasTkAgg(fig, master=self.strategy_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas.draw()

    def create_track_page(self):
        ctk.CTkLabel(self.content_area, text="Track Visualization", font=("Roboto", 24, "bold")).pack(pady=10)
        
        btn_frame = ctk.CTkFrame(self.content_area)
        btn_frame.pack(pady=10)
        
        self.track_session_var = ctk.StringVar(value="Race")
        ctk.CTkOptionMenu(btn_frame, variable=self.track_session_var, values=["Race", "Qualifying"], width=100, fg_color="#333", button_color="#444").pack(side="left", padx=5)
        
        ctk.CTkButton(btn_frame, text="Gear Shift Map", command=self.plot_gear_map, fg_color="#E10600").pack(side="left", padx=10)
        ctk.CTkButton(btn_frame, text="Speed Map", command=self.plot_speed_map, fg_color="#2CC985").pack(side="left", padx=10)
        
        # Track Dominance Inputs
        dom_frame = ctk.CTkFrame(self.content_area)
        dom_frame.pack(pady=5)
        
        self.dom_d1_var = ctk.StringVar(value="VER")
        self.dom_d2_var = ctk.StringVar(value="HAM")
        
        ctk.CTkLabel(dom_frame, text="Dominance:").pack(side="left", padx=5)
        ctk.CTkEntry(dom_frame, textvariable=self.dom_d1_var, width=50).pack(side="left", padx=5)
        ctk.CTkLabel(dom_frame, text="vs").pack(side="left", padx=5)
        ctk.CTkEntry(dom_frame, textvariable=self.dom_d2_var, width=50).pack(side="left", padx=5)
        
        ctk.CTkButton(dom_frame, text="Show Dominance", command=self.plot_track_dominance, fg_color="#FF8700").pack(side="left", padx=10)
        ctk.CTkButton(dom_frame, text="Battle Map (Gap)", command=self.plot_battle_map, fg_color="#9B59B6").pack(side="left", padx=10)
        
        self.track_frame = ctk.CTkFrame(self.content_area, fg_color="transparent")
        self.track_frame.pack(fill="both", expand=True, padx=10, pady=10)

    def plot_gear_map(self):
        session_type = self.track_session_var.get()
        session = self.load_session(session_type)
        if not session: return
        
        # Pick fastest lap of the session
        lap = session.laps.pick_fastest()
        tel = lap.get_telemetry()
        
        if 'nGear' not in tel.columns: return
        
        x = np.array(tel['X'].values)
        y = np.array(tel['Y'].values)
        
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        gear = tel['nGear'].to_numpy().astype(float)
        
        import matplotlib.collections as mcoll
        
        for widget in self.track_frame.winfo_children(): widget.destroy()
        fig, ax = plt.subplots(figsize=(8, 6), facecolor='#121212')
        ax.set_facecolor('#121212')
        
        cmap = plt.get_cmap('Paired', 8)
        lc = mcoll.LineCollection(segments, cmap=cmap, norm=plt.Normalize(1, 8))
        lc.set_array(gear)
        lc.set_linewidth(4)
        
        ax.add_collection(lc)
        ax.axis('off')
        ax.set_aspect('equal')
        
        cbar = plt.colorbar(lc, ax=ax, ticks=np.arange(1, 9))
        cbar.ax.tick_params(colors='white')
        cbar.set_label('Gear', color='white')
        
        ax.set_title(f"Gear Shift Map - {lap['Driver']}", color="white")
        ax.autoscale_view()
        
        canvas = FigureCanvasTkAgg(fig, master=self.track_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas.draw()

    def plot_speed_map(self):
        # Similar but Speed
        session_type = self.track_session_var.get()
        session = self.load_session(session_type)
        if not session: return
        lap = session.laps.pick_fastest()
        tel = lap.get_telemetry()
        
        x = np.array(tel['X'].values)
        y = np.array(tel['Y'].values)
        speed = tel['Speed'].to_numpy()
        
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        import matplotlib.collections as mcoll
        
        for widget in self.track_frame.winfo_children(): widget.destroy()
        fig, ax = plt.subplots(figsize=(8, 6), facecolor='#121212')
        ax.set_facecolor('#121212')
        
        lc = mcoll.LineCollection(segments, cmap='plasma', norm=plt.Normalize(speed.min(), speed.max()))
        lc.set_array(speed)
        lc.set_linewidth(4)
        ax.add_collection(lc)
        ax.axis('off')
        ax.set_aspect('equal')
        
        cbar = plt.colorbar(lc, ax=ax)
        cbar.ax.tick_params(colors='white')
        cbar.set_label('Speed (km/h)', color='white')
        
        ax.set_title(f"Track Speed Map - {lap['Driver']}", color="white")
        ax.autoscale_view()
        
        canvas = FigureCanvasTkAgg(fig, master=self.track_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas.draw()

    def plot_track_dominance(self):
        session_type = self.track_session_var.get()
        session = self.load_session(session_type)
        if not session: return
        
        d1 = self.dom_d1_var.get().upper()
        d2 = self.dom_d2_var.get().upper()
        
        try:
            laps_d1 = session.laps.pick_drivers(d1).pick_fastest()
            laps_d2 = session.laps.pick_drivers(d2).pick_fastest()
            
            if laps_d1 is None or laps_d2 is None:
                messagebox.showerror("Error", "Brief data not available for one or both drivers.")
                return

            # Telemetry
            tel_d1 = laps_d1.get_telemetry().add_distance()
            tel_d2 = laps_d2.get_telemetry().add_distance()
            
            # Create a common distance array
            total_dist = max(tel_d1['Distance'].max(), tel_d2['Distance'].max())
            dist = np.linspace(0, total_dist, num=100) # 100 minisectors
            
            # Interpolate Time
            t1 = np.interp(dist, tel_d1['Distance'], tel_d1['Time'].dt.total_seconds())
            t2 = np.interp(dist, tel_d2['Distance'], tel_d2['Time'].dt.total_seconds())
            
            # Calculate delta for each segment
            dt1 = np.diff(t1)
            dt2 = np.diff(t2)
            
            # Determine dominant driver for each segment
            dominance = []
            for i in range(len(dt1)):
                if dt1[i] < dt2[i]:
                    dominance.append(1) # D1
                else:
                    dominance.append(2) # D2
                    
            # Get coords for segments
            x_track = tel_d1['X'].to_numpy()
            y_track = tel_d1['Y'].to_numpy()
            dist_track = tel_d1['Distance'].to_numpy()
            
            x_interp = np.interp(dist, dist_track, x_track)
            y_interp = np.interp(dist, dist_track, y_track)
            
            points = np.array([x_interp, y_interp]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            # Plot
            for widget in self.track_frame.winfo_children(): widget.destroy()
            fig, ax = plt.subplots(figsize=(8, 6), facecolor='#121212')
            ax.set_facecolor('#121212')
            
            import matplotlib.collections as mcoll
            from matplotlib.colors import ListedColormap, BoundaryNorm
            
            c1 = fastf1.plotting.get_driver_color(d1, session=session)
            c2 = fastf1.plotting.get_driver_color(d2, session=session)
            
            # Create custom colormap
            cmap = ListedColormap([c1, c2])
            norm = BoundaryNorm([0.5, 1.5, 2.5], cmap.N)
            
            lc = mcoll.LineCollection(segments, cmap=cmap, norm=norm)
            lc.set_array(np.array(dominance))
            lc.set_linewidth(5)
            
            ax.add_collection(lc)
            ax.axis('off')
            ax.set_aspect('equal')
            
            # Custom Legend
            legend_lines = [
                plt.Line2D([0], [0], color=c1, linewidth=4, label=d1),
                plt.Line2D([0], [0], color=c2, linewidth=4, label=d2)
            ]
            ax.legend(handles=legend_lines, facecolor='#333', labelcolor='white')
            
            ax.set_title(f"Track Dominance: {d1} vs {d2} ({session_type})", color="white")
            ax.autoscale_view()
            
            canvas = FigureCanvasTkAgg(fig, master=self.track_frame)
            canvas.get_tk_widget().pack(fill="both", expand=True)
            canvas.draw()
            
        except Exception as e:
            messagebox.showerror("Error", f"Analysis Failed: {e}")
            print(e)
            return

    def plot_battle_map(self):
        """
        Visualizes where drivers gained/lost time against each other (Gap Derivative).
        """
        session_type = self.track_session_var.get()
        session = self.load_session(session_type)
        if not session: return
        
        d1 = self.dom_d1_var.get().upper()
        d2 = self.dom_d2_var.get().upper()
        
        try:
            laps_d1 = session.laps.pick_drivers(d1).pick_fastest()
            laps_d2 = session.laps.pick_drivers(d2).pick_fastest()
            
            if laps_d1 is None or laps_d2 is None:
                messagebox.showerror("Error", "Data not available.")
                return

            tel_d1 = laps_d1.get_telemetry().add_distance()
            tel_d2 = laps_d2.get_telemetry().add_distance()
            
            # 1. Align Data
            total_dist = max(tel_d1['Distance'].max(), tel_d2['Distance'].max())
            dist = np.linspace(0, total_dist, num=400) # Higher res for smooth gradient
            
            t1 = np.interp(dist, tel_d1['Distance'], tel_d1['Time'].dt.total_seconds())
            t2 = np.interp(dist, tel_d2['Distance'], tel_d2['Time'].dt.total_seconds())
            
            # 2. Calculate Gap and Gradient (Time Gained per Distance)
            gap = t2 - t1 # +ve if D1 is ahead (t2 > t1)
            
            # We want to show "Change in Gap" -> Who is faster HERE?
            # derivative of gap w.r.t distance
            gap_deriv = np.gradient(gap)
            
            # Smoothing
            from scipy.ndimage import gaussian_filter1d
            gap_deriv = gaussian_filter1d(gap_deriv, sigma=5)
            
            # Normalize for color (symmetric around 0)
            # +ve deriv: Gap increasing (D1 is faster)
            # -ve deriv: Gap decreasing (D2 is faster)
            limit = max(abs(gap_deriv.min()), abs(gap_deriv.max()))
            
            # 3. Track coords
            x_track = tel_d1['X'].to_numpy()
            y_track = tel_d1['Y'].to_numpy()
            dist_track = tel_d1['Distance'].to_numpy()
            
            x_interp = np.interp(dist, dist_track, x_track)
            y_interp = np.interp(dist, dist_track, y_track)
            
            points = np.array([x_interp, y_interp]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            
            # 4. Plot
            for widget in self.track_frame.winfo_children(): widget.destroy()
            fig, ax = plt.subplots(figsize=(8, 6), facecolor='#121212')
            ax.set_facecolor('#121212')
            
            import matplotlib.collections as mcoll
            from matplotlib.colors import LinearSegmentedColormap
            
            c1 = fastf1.plotting.get_driver_color(d1, session=session)
            c2 = fastf1.plotting.get_driver_color(d2, session=session)
            
            # Custom divergent colormap: D2 (Faster) -> D1 (Faster)
            # D2 is negative deriv, D1 is positive deriv
            colors = [c2, '#222222', c1] # Blue -> Grey -> Red
            cmap = LinearSegmentedColormap.from_list("BattleMap", colors)
            
            # Norm centered at 0
            norm = plt.Normalize(-limit, limit)
            
            lc = mcoll.LineCollection(segments, cmap=cmap, norm=norm)
            lc.set_array(gap_deriv)
            lc.set_linewidth(6)
            ax.add_collection(lc)
            
            ax.axis('off')
            ax.set_aspect('equal')
            
            # Legend
            legend_lines = [
                plt.Line2D([0], [0], color=c1, linewidth=4, label=f"{d1} Faster"),
                plt.Line2D([0], [0], color=c2, linewidth=4, label=f"{d2} Faster")
            ]
            ax.legend(handles=legend_lines, facecolor='#333', labelcolor='white')
            
            ax.set_title(f"Battle Map: where time was won/lost\n{d1} vs {d2}", color="white")
            ax.autoscale_view()
            
            canvas = FigureCanvasTkAgg(fig, master=self.track_frame)
            canvas.get_tk_widget().pack(fill="both", expand=True)
            canvas.draw()
            
        except Exception as e:
            messagebox.showerror("Error", f"Analysis Failed: {e}")

    def create_dna_page(self):
        ctk.CTkLabel(self.content_area, text="Driver DNA", font=("Roboto", 24, "bold")).pack(pady=10)
        ctk.CTkButton(self.content_area, text="Generate DNA", command=self.plot_dna, fg_color="#E10600").pack(pady=5)
        self.dna_frame = ctk.CTkFrame(self.content_area, fg_color="transparent")
        self.dna_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
    def plot_dna(self):
        session = self.load_session()
        if not session: return
        
        drivers = ['VER', 'HAM', 'NOR', 'LEC', 'ALO'] # Compare top 5 for clarity
        
        stats = {'Driver': [], 'Quali': [], 'RacePace': [], 'Consistency': [], 'TyreMgmt': []}
        
        pole = session.laps.pick_fastest()['LapTime'].total_seconds()
        
        for d in drivers:
            laps = session.laps.pick_drivers(d)
            if laps.empty: continue
            
            # Quali (Gap to Pole inverted? Or just raw gap normalized?)
            # Let's map Gap [0, 2s] to [100, 0]
            best_lap = laps.pick_fastest()['LapTime'].total_seconds()
            gap = best_lap - pole
            quali_score = max(0, 100 - (gap * 50))
            
            # Race Pace (Avg of quick laps)
            quick_laps = laps.pick_quicklaps()
            if quick_laps.empty: continue
            avg_pace = quick_laps['LapTime'].dt.total_seconds().mean()
            # Normalize: Pole vs Avg Pace
            pace_gap = avg_pace - pole
            pace_score = max(0, 100 - (pace_gap * 20))
            
            # Consistency (1 / StdDev)
            std = quick_laps['LapTime'].dt.total_seconds().std()
            cons_score = max(0, 100 - (std * 50))
            
            # Tyre Mgmt (Avg Stint Length)
            stints = laps.groupby('Stint')['LapNumber'].count().mean()
            tyre_score = min(100, stints * 5) # 20 laps avg = 100
            
            stats['Driver'].append(d)
            stats['Quali'].append(quali_score)
            stats['RacePace'].append(pace_score)
            stats['Consistency'].append(cons_score)
            stats['TyreMgmt'].append(tyre_score)
            
        # Radar Chart
        categories = ['Quali', 'RacePace', 'Consistency', 'TyreMgmt']
        N = len(categories)
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]
        
        for widget in self.dna_frame.winfo_children(): widget.destroy()
        fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True), facecolor='#121212')
        ax.set_facecolor('#121212')
        
        for i, d in enumerate(stats['Driver']):
            values = [stats[c][i] for c in categories]
            values += values[:1]
            c = fastf1.plotting.get_driver_color(d, session=session)
            ax.plot(angles, values, linewidth=2, linestyle='solid', label=d, color=c)
            
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, color="white")
        ax.tick_params(colors="white")
        ax.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1), facecolor='#333', labelcolor='white')
        ax.set_title("Driver DNA Comparison", color="white", y=1.1)
        
        canvas = FigureCanvasTkAgg(fig, master=self.dna_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas.draw()
        
    def create_position_page(self):
        ctk.CTkLabel(self.content_area, text="Position Change Chart", font=("Roboto", 24, "bold")).pack(pady=10)
        ctk.CTkButton(self.content_area, text="Load Chart", command=self.plot_position_chart, fg_color="#E10600").pack(pady=5)
        self.pos_frame = ctk.CTkFrame(self.content_area, fg_color="transparent")
        self.pos_frame.pack(fill="both", expand=True, padx=10, pady=10)

    def plot_position_chart(self):
        session = self.load_session()
        if not session: return
        
        for widget in self.pos_frame.winfo_children(): widget.destroy()
        fig, ax = plt.subplots(figsize=(10, 6), facecolor='#121212')
        ax.set_facecolor('#121212')
        
        drivers = session.drivers
        # Fallback colors if FastF1 returns white/none or fails
        fallback_colors = plt.cm.tab20.colors 
        
        for i, drv in enumerate(drivers):
             laps = session.laps.pick_drivers(drv)
             try:
                 color = fastf1.plotting.get_driver_color(drv, session=session)
             except:
                 color = None
            
             # If color is invalid or white (often default fallback), use our fallback list
             if not color or color.lower() in ['white', '#ffffff']:
                 color = fallback_colors[i % len(fallback_colors)]
                 
             ax.plot(laps['LapNumber'], laps['Position'], color=color, label=drv)
             
        ax.set_ylabel("Position", color="white")
        ax.set_xlabel("Lap", color="white")
        ax.invert_yaxis()
        ax.tick_params(colors="white")
        ax.set_title("Position Changes", color="white")
        # Add legend outside the plot to avoid clutter
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., 
                  facecolor='#333', labelcolor='white', fontsize='small')
        fig.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig, master=self.pos_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas.draw()

    def create_race_lines_page(self):
        ctk.CTkLabel(self.content_area, text="Race Line Comparison", font=("Roboto", 24, "bold")).pack(pady=10)
        
        controls = ctk.CTkFrame(self.content_area)
        controls.pack(pady=10)
        
        self.rl_d1_var = ctk.StringVar(value="VER")
        self.rl_d2_var = ctk.StringVar(value="HAM")
        
        ctk.CTkEntry(controls, textvariable=self.rl_d1_var, width=60).pack(side="left", padx=10)
        ctk.CTkLabel(controls, text="VS").pack(side="left", padx=10)
        ctk.CTkEntry(controls, textvariable=self.rl_d2_var, width=60).pack(side="left", padx=10)
        
        ctk.CTkButton(controls, text="Analyze Race Lines", command=self.plot_race_lines, fg_color="#E10600").pack(side="left", padx=20)
        
        ctk.CTkLabel(self.content_area, text="Click on track to zoom in/out", font=("Roboto", 12), text_color="gray").pack()

        self.race_lines_frame = ctk.CTkFrame(self.content_area, fg_color="transparent")
        self.race_lines_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
    def smooth_track_data(self, x, y, s=None, n_points=None):
        """Smooths track data using B-spline interpolation."""
        if len(x) < 4: return x, y
        if n_points is None: n_points = len(x)
        try:
             # Remove duplicates for splprep
             points = np.vstack((x, y)).T
             _, idx = np.unique(points, axis=0, return_index=True)
             idx.sort()
             x = x[idx]
             y = y[idx]
             if len(x) < 4: return x, y
             
             tck, u = splprep([x, y], s=s, per=0) 
             u_new = np.linspace(u.min(), u.max(), n_points)
             x_new, y_new = splev(u_new, tck)
             return x_new, y_new
        except Exception as e:
             print(f"Smoothing failed: {e}")
             return x, y

    def generate_track_ribbon(self, x, y, width=120):
        """Generates inner and outer track edges with curvature-based shifting."""
        # Calculate gradients
        dx = np.gradient(x)
        dy = np.gradient(y)
        ddx = np.gradient(dx)
        ddy = np.gradient(dy)
        
        # Normals
        ds = np.sqrt(dx**2 + dy**2)
        ds[ds==0] = 1.0 
        
        nx = -dy / ds
        ny = dx / ds
        
        # Calculate Curvature (signed)
        # k = (dx*ddy - dy*ddx) / (dx^2 + dy^2)^1.5
        k = (dx * ddy - dy * ddx) / (ds**3)
        
        # Shift Logic:
        # If k > 0 (Left Turn), center of curvature is Left (Normal direction).
        # We want the Racing Line (x,y) to be on the Left edge (Apex).
        # So we need to shift the Track Center to the Right (Opposite Normal).
        # So shift should be negative proportional to k.
        
        # Limit shift to slightly less than half width to ensure line stays within track
        max_shift = width * 0.45 
        
        # Scale factor: tuned manually. 
        # k values are small (e.g. 1e-3). Width is ~80.
        # We want k~high to give shift~max_shift.
        # Let's normalize k roughly.
        k_max_est = np.max(np.abs(k))
        if k_max_est == 0: k_max_est = 1.0
        
        shift_factor = k / k_max_est * max_shift * -1.2 # Boost a bit to ensure hitting apex
        
        # Apply clamping
        shift_factor = np.clip(shift_factor, -max_shift, max_shift)
        
        # Shifted Center
        x_c = x + nx * shift_factor
        y_c = y + ny * shift_factor
        
        # Edges from new center
        x_out = x_c + nx * width * 0.5
        y_out = y_c + ny * width * 0.5
        
        x_in = x_c - nx * width * 0.5
        y_in = y_c - ny * width * 0.5
        
        return x_out, y_out, x_in, y_in

    def plot_race_lines(self):
        session = self.load_session()
        if not session: return
        
        d1 = self.rl_d1_var.get().upper()
        d2 = self.rl_d2_var.get().upper()
        
        try:
            laps_d1 = session.laps.pick_drivers(d1).pick_fastest()
            laps_d2 = session.laps.pick_drivers(d2).pick_fastest()
            
            tel_d1 = laps_d1.get_telemetry()
            tel_d2 = laps_d2.get_telemetry()
            
            # Get track reference
            fastest_lap = session.laps.pick_fastest()
            tel_track = fastest_lap.get_telemetry()
            
            x_track, y_track = tel_track['X'].to_numpy(), tel_track['Y'].to_numpy()
            
            # 1. Smooth the track center line
            x_track_smooth, y_track_smooth = self.smooth_track_data(x_track, y_track, s=10000, n_points=len(x_track)*5)
            
            # 2. Generate Ribbon with Apex Logic
            x_out, y_out, x_in, y_in = self.generate_track_ribbon(x_track_smooth, y_track_smooth, width=90)
            
            # 3. Smooth Driver Lines
            x_d1, y_d1 = tel_d1['X'].to_numpy(), tel_d1['Y'].to_numpy()
            x_d2, y_d2 = tel_d2['X'].to_numpy(), tel_d2['Y'].to_numpy()
            
            x_d1_s, y_d1_s = self.smooth_track_data(x_d1, y_d1, s=0, n_points=len(x_d1)*3)
            x_d2_s, y_d2_s = self.smooth_track_data(x_d2, y_d2, s=0, n_points=len(x_d2)*3)
            
        except Exception as e:
            messagebox.showerror("Error", f"Could not load data for drivers: {e}")
            return

        c1 = fastf1.plotting.get_driver_color(d1, session=session)
        c2 = fastf1.plotting.get_driver_color(d2, session=session)
        
        for widget in self.race_lines_frame.winfo_children(): widget.destroy()
        
        fig, ax = plt.subplots(figsize=(10, 8), facecolor='#121212')
        ax.set_facecolor('#121212')
        
        # Plot Track Ribbon
        poly_x = np.concatenate([x_out, x_in[::-1]])
        poly_y = np.concatenate([y_out, y_in[::-1]])
        
        # Use a nice dark gray for asphalt
        ax.fill(poly_x, poly_y, color='#252525', zorder=1, label="Track Surface")
        
        # Optional: Plot borders
        ax.plot(x_out, y_out, color='#444444', linewidth=0.5, zorder=1.5)
        ax.plot(x_in, y_in, color='#444444', linewidth=0.5, zorder=1.5)
        
        # Plot Driver Lines
        ax.plot(x_d1_s, y_d1_s, color=c1, label=d1, linewidth=2, zorder=2)
        ax.plot(x_d2_s, y_d2_s, color=c2, label=d2, linewidth=2, linestyle='--', zorder=3)
        
        ax.set_title(f"Race Line Comparison: {d1} vs {d2}", color="white")
        ax.axis('equal')
        ax.axis('off')
        
        # Custom Legend
        legend_lines = [
            plt.Line2D([0], [0], color=c1, linewidth=2, label=d1),
            plt.Line2D([0], [0], color=c2, linewidth=2, linestyle='--', label=d2),
            plt.Rectangle((0,0), 1, 1, facecolor='#252525', edgecolor='#444444', label='Track')
        ]
        ax.legend(handles=legend_lines, facecolor='#333', labelcolor='white')
        
        self.rl_ax = ax
        self.rl_is_zoomed = False
        self.rl_fig = fig
        
        canvas = FigureCanvasTkAgg(fig, master=self.race_lines_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas.draw()
        
        self.rl_canvas = canvas
        
        # Click event for zoom
        def on_click(event):
            if event.inaxes != self.rl_ax: return
            if event.button != 1: return # Left click only
            
            if self.rl_is_zoomed:
                self.rl_ax.autoscale()
                self.rl_is_zoomed = False
            else:
                x, y = event.xdata, event.ydata
                # Zoom radius
                zoom_radius = 600 
                self.rl_ax.set_xlim(x - zoom_radius, x + zoom_radius)
                self.rl_ax.set_ylim(y - zoom_radius, y + zoom_radius)
                self.rl_is_zoomed = True
                
            self.rl_canvas.draw()
            
        fig.canvas.mpl_connect('button_press_event', on_click)

    def create_tyre_page(self):
        ctk.CTkLabel(self.content_area, text="Tyre Degradation Analysis", font=("Roboto", 24, "bold")).pack(pady=10)
        
        controls = ctk.CTkFrame(self.content_area)
        controls.pack(pady=10)
        
        self.tyre_d1_var = ctk.StringVar(value="VER")
        self.tyre_d2_var = ctk.StringVar(value="HAM")
        
        ctk.CTkEntry(controls, textvariable=self.tyre_d1_var, width=60).pack(side="left", padx=10)
        ctk.CTkLabel(controls, text="VS").pack(side="left", padx=10)
        ctk.CTkEntry(controls, textvariable=self.tyre_d2_var, width=60).pack(side="left", padx=10)
        
        ctk.CTkButton(controls, text="Analyze Degradation", command=self.plot_tyre_degradation, fg_color="#E10600").pack(side="left", padx=20)
        
        self.tyre_frame = ctk.CTkFrame(self.content_area, fg_color="transparent")
        self.tyre_frame.pack(fill="both", expand=True, padx=10, pady=10)

    def plot_tyre_degradation(self):
        session = self.load_session()
        if not session: return
        
        d1 = self.tyre_d1_var.get().upper()
        d2 = self.tyre_d2_var.get().upper()
        
        try:
            # Get laps for both drivers
            laps_d1 = session.laps.pick_drivers(d1)
            laps_d2 = session.laps.pick_drivers(d2)
            
            # Filter for quicklaps to remove in-laps/out-laps/SC
            laps_d1 = laps_d1.pick_quicklaps(1.07)
            laps_d2 = laps_d2.pick_quicklaps(1.07)
            
            # Ensure we have TyreLife
            if 'TyreLife' not in laps_d1.columns:
                pass 
                
        except Exception as e:
            messagebox.showerror("Error", f"Data Error: {e}")
            return
            
        for widget in self.tyre_frame.winfo_children(): widget.destroy()
        
        fig, ax = plt.subplots(figsize=(10, 8), facecolor='#121212')
        ax.set_facecolor('#121212')
        
        # Plot D1
        for stint in laps_d1['Stint'].unique():
            stint_laps = laps_d1[laps_d1['Stint'] == stint]
            if stint_laps.empty: continue
            
            compound = stint_laps['Compound'].iloc[0]
            color = fastf1.plotting.get_driver_color(d1, session=session)
            comp_color = self.get_compound_color_safe(compound, session=session)
            
            ax.scatter(stint_laps['TyreLife'], stint_laps['LapTime'].dt.total_seconds(), 
                       color=comp_color, marker='o', edgecolors=color, linewidth=2, s=50, label=f"{d1} {compound}")

        # Plot D2
        for stint in laps_d2['Stint'].unique():
            stint_laps = laps_d2[laps_d2['Stint'] == stint]
            if stint_laps.empty: continue
            
            compound = stint_laps['Compound'].iloc[0]
            color = fastf1.plotting.get_driver_color(d2, session=session)
            comp_color = self.get_compound_color_safe(compound, session=session)
            
            ax.scatter(stint_laps['TyreLife'], stint_laps['LapTime'].dt.total_seconds(), 
                       color=comp_color, marker='x', linewidth=2, s=50, label=f"{d2} {compound}")

        ax.set_xlabel("Tyre Age (Laps)", color="white")
        ax.set_ylabel("Lap Time (s)", color="white")
        ax.tick_params(colors="white")
        ax.set_title(f"Tyre Degradation: {d1} (o) vs {d2} (x)", color="white")
        
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), facecolor='#333', labelcolor='white')
        
        canvas = FigureCanvasTkAgg(fig, master=self.tyre_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas.draw()

    def add_placeholder_plot(self, title):
        fig, ax = plt.subplots(figsize=(8, 5), facecolor='#121212')
        ax.set_facecolor('#121212')
        ax.text(0.5, 0.5, f"{title}\n\nSelect Drivers & Session to Start Analysis", 
                color='gray', ha='center', va='center', fontsize=14, fontfamily="Roboto")
        ax.axis('off')
        
        canvas = FigureCanvasTkAgg(fig, master=self.content_area)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=20, pady=20)
        canvas.draw()

    def on_closing(self):
        self.is_running = False
        if hasattr(self, 'anim_running'):
            self.anim_running = False
        self.quit()
        self.destroy()
        sys.exit()
    # --- 3D CAMERA INTERACTION METHODS ---
    def on_cam_press(self, event):
        if not hasattr(self, 'track_3d_var') or not self.track_3d_var.get(): return
        self.is_dragging = True
        self.last_mouse_x = event.x
        self.last_mouse_y = event.y

    def on_cam_drag(self, event):
        if not hasattr(self, 'track_3d_var') or not self.track_3d_var.get() or not self.is_dragging: return
        
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
        # Only run if window is open
        if not self.is_running: return

        if hasattr(self, 'track_3d_var') and self.track_3d_var.get() and hasattr(self, 'track_ax') and self.track_ax:
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
                     self.track_ax.view_init(elev=self.cam_elev, azim=self.cam_azim)
                     self.track_canvas.draw_idle()
                 except: pass

        # Run loop at ~60FPS (16ms)
        self.after(16, self.update_camera_physics)

    def on_scroll(self, event):
        """ Zoom in/out with scroll wheel in 3D mode by adjusting axis limits """
        if not hasattr(self, 'track_3d_var') or not self.track_3d_var.get():
            return
            
        if not hasattr(self, 'track_ax') or event.inaxes != self.track_ax:
            return

        # Zoom Factor logic
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
        self.track_canvas.draw_idle()

    def apply_zoom(self):
        """ Apply current zoom level to axis limits """
        if not hasattr(self, 'plot_bounds') or not hasattr(self, 'track_ax'):
            return
            
        (min_x, max_x, min_y, max_y, min_z, max_z) = self.plot_bounds
        
        center_x = (min_x + max_x) / 2
        center_y = (min_y + max_y) / 2
        
        span_x = (max_x - min_x)
        span_y = (max_y - min_y)
        
        # Apply Zoom
        new_span_x = span_x * self.zoom_level
        new_span_y = span_y * self.zoom_level
        
        self.track_ax.set_xlim(center_x - new_span_x/2, center_x + new_span_x/2)
        self.track_ax.set_ylim(center_y - new_span_y/2, center_y + new_span_y/2)

    def create_track_info_page(self):
        ctk.CTkLabel(self.content_area, text="Detailed Track Info", font=("Roboto", 24, "bold")).pack(pady=10)
        
        controls = ctk.CTkFrame(self.content_area)
        controls.pack(pady=10)
        
        self.track_session_var = ctk.StringVar(value="Race")
        ctk.CTkOptionMenu(controls, variable=self.track_session_var, values=["Race", "Qualifying"], width=100, fg_color="#333", button_color="#444").pack(side="left", padx=5)
        
        ctk.CTkButton(controls, text="Generate Map", command=self.plot_track_info, fg_color="#3498DB").pack(side="left", padx=10)
        
        # 3D Toggle
        self.track_3d_var = ctk.BooleanVar(value=True) # Default to 3D for impact
        self.switch_3d = ctk.CTkSwitch(controls, text="3D View", variable=self.track_3d_var, onvalue=True, offvalue=False, command=self.plot_track_info)
        self.switch_3d.pack(side="left", padx=20)
        
        self.track_frame = ctk.CTkFrame(self.content_area, fg_color="transparent")
        self.track_frame.pack(fill="both", expand=True, padx=10, pady=10)
        # REMOVED 3D Toggle Switch
        
        self.track_frame = ctk.CTkFrame(self.content_area, fg_color="#121212")
        self.track_frame.pack(fill="both", expand=True, padx=20, pady=10)
        
        self.add_placeholder_plot("Track Map Area")

    def plot_track_info(self):
        """ Plot High-Quality 2D Track Map (Double Line Style) """
        session_type = self.track_session_var.get()
        session = self.load_session(session_type)
        if not session: return
        
        try:
            lap = session.laps.pick_fastest()
            tel = lap.get_telemetry()
            
            x = np.array(tel['X'].values)
            y = np.array(tel['Y'].values)
            drs = tel['DRS'].to_numpy() # 0=Off, 1=?, 8+=On
            
            # 2. Sector Calculations
            start_time = lap['LapStartTime']
            s1_end_time = start_time + lap['Sector1Time']
            s2_end_time = s1_end_time + lap['Sector2Time']
            
            times = tel['Time']
            idx_s1 = tel.index.get_loc((times - s1_end_time).abs().idxmin())
            idx_s2 = tel.index.get_loc((times - s2_end_time).abs().idxmin())
            
            # Clear previous
            for widget in self.track_frame.winfo_children(): widget.destroy()
            
            # Zone Name
            year_val = 2023
            try: year_val = int(self.year_var.get())
            except: pass
            zone_label = "Overtake Zone" if year_val >= 2026 else "DRS Zone"
            
            # --- 2D PLOT SETUP ---
            fig, ax = plt.subplots(figsize=(10, 7), facecolor='black')
            ax.set_facecolor('black')
            
            # Remove Margins to maximize map size
            plt.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.02)
            
            # --- TRACK STYLE: DOUBLE LINE ---
            # 1. Outer Border (Grey)
            ax.plot(x, y, color='#AAAAAA', linewidth=10, zorder=1)
            
            # 2. Inner Fill (Black)
            ax.plot(x, y, color='black', linewidth=6, zorder=2)
            
            # --- ZONES (DRS/Overtake) ---
            # Green "In-fill" ON TOP of black fill
            drs_active = drs > 8
            if drs_active.sum() > 0:
                drs_x = x.copy(); drs_y = y.copy()
                drs_x[~drs_active] = np.nan
                drs_y[~drs_active] = np.nan
                
                # Plot GREEN on top of black fill (zorder=3)
                # Width 5 fits inside 6, but let's match it (6) or go slightly thinner (4)
                # If zorder=3, it overlays black.
                ax.plot(drs_x, drs_y, color='#00FF00', linewidth=5, label=zone_label, zorder=3)
            
            # --- SECTOR MARKINGS (Clean with Offset) ---
            def draw_sector_marker(idx, label, color='white', offset=(20, 20)):
                lx, ly = x[idx], y[idx]
                
                # 1. MARKER on track
                ax.scatter([lx], [ly], color='white', s=60, marker='s', zorder=4, edgecolor='black', linewidth=1)
                
                # 2. LABEL with Arrow pointing to marker
                # Use annotate to place text away from track
                ax.annotate(label, 
                            xy=(lx, ly), 
                            xytext=offset, 
                            textcoords='offset points',
                            color=color, fontsize=10, fontweight='bold',
                            ha='center', va='center', zorder=5,
                            arrowprops=dict(facecolor=color, arrowstyle="-", linewidth=1.5),
                            bbox=dict(facecolor='black', edgecolor=color, boxstyle='round,pad=0.2', alpha=0.8))

            # Offsets can be dynamic, but fixed is safer than overlapping track.
            # Try to push them "outwards".
            # For now, distinct offsets.
            draw_sector_marker(0, "FINISH", color='#FFFFFF', offset=(30, 0))
            draw_sector_marker(idx_s1, "S1", color='#FFD700', offset=(-30, 30))
            draw_sector_marker(idx_s2, "S2", color='#FFD700', offset=(30, -30))
            
            # --- CORNER NUMBERS ---
            # Keep them on track or slightly off? On track is standard.
            if hasattr(session, 'circuit_info'):
                try:
                    corners = session.circuit_info.corners
                    for _, corner in corners.iterrows():
                        ax.text(corner['X'], corner['Y'], f"{corner['Number']}", 
                                color='#888888', fontsize=7, ha='center', va='center', fontweight='bold', zorder=2.5)
                except: pass

            # Layout
            ax.axis('off')
            ax.set_aspect('equal')
            ax.set_title(f"{session.event.EventName} - Track Map ({year_val})", color='white', fontsize=12)
            
            # Legend
            if drs_active.sum() > 0:
                leg = ax.legend(facecolor='black', edgecolor='#444', loc='upper right', fontsize=8)
                for text in leg.get_texts(): text.set_color("white")
            
            # Embed
            canvas = FigureCanvasTkAgg(fig, master=self.track_frame)
            self.track_canvas = canvas
            canvas.get_tk_widget().pack(fill="both", expand=True)
            canvas.draw()
            
        except Exception as e:
            print(f"Track Info Error: {e}")
            import traceback
            traceback.print_exc()

    def create_track_info_page(self):
        # Header
        ctk.CTkLabel(self.content_area, text="Detailed Track Info", font=("Arial", 24, "bold")).pack(pady=10)
        
        controls = ctk.CTkFrame(self.content_area)
        controls.pack(pady=10)
        
        self.track_session_var = ctk.StringVar(value="Race")
        ctk.CTkOptionMenu(controls, variable=self.track_session_var, values=["Race", "Qualifying"], width=100, fg_color="#333", button_color="#444").pack(side="left", padx=5)
        
        self.year_var = ctk.StringVar(value="2023")  
        
        ctk.CTkButton(controls, text="Generate Map", command=self.plot_track_info, fg_color="#3498DB").pack(side="left", padx=10)
        
        self.track_frame = ctk.CTkFrame(self.content_area, fg_color="#121212")
        self.track_frame.pack(fill="both", expand=True, padx=20, pady=10)
        
        # Placeholder
        fig, ax = plt.subplots(figsize=(8, 5), facecolor='#121212')
        ax.set_facecolor('#121212')
        ax.text(0.5, 0.5, "Select Session & Click Generate", color='gray', ha='center', va='center', fontsize=14)
        ax.axis('off')
        
        canvas = FigureCanvasTkAgg(fig, master=self.track_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=20, pady=20)
        canvas.draw()

if __name__ == "__main__":
    import sys
    app = AnalyticsDashboardApp()
    app.mainloop()