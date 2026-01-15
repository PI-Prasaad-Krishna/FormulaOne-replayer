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

# Configure appearance
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

class AnalyticsDashboardApp(ctk.CTkToplevel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.title("F1 Analytics Dashboard 2025")
        self.geometry("1400x900")
        
        # Grid Layout
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

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

        # --- Main Content Area ---
        self.content_area = ctk.CTkFrame(self, corner_radius=0, fg_color="#121212")
        self.content_area.grid(row=0, column=1, sticky="nsew")
        
        # Default Page
        self.show_page("telemetry")

    def show_page(self, page_name):
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

    def load_session(self):
        # Helper to load session if not loaded (hardcoded 2023 Abu Dhabi for demo or add inputs)
        # For simplicity, we'll use a dialog or just 2023 Abu Dhabi
        try:
             year = int(self.year_var.get())
             circuit = self.event_var.get()
        except:
             messagebox.showerror("Error", "Invalid Year")
             return None
             
        status_lbl = ctk.CTkLabel(self.sidebar, text="Loading Data...", text_color="orange")
        status_lbl.pack(side="bottom", pady=10)
        self.update()
        
        try:
            # Enable Cache
            if not os.path.exists("f1_cache"): os.makedirs("f1_cache")
            fastf1.Cache.enable_cache("f1_cache")
            
            session = fastf1.get_session(year, circuit, 'R')
            session.load(telemetry=True, laps=True, weather=False)
            status_lbl.destroy()
            return session
        except Exception as e:
             status_lbl.configure(text="Error", text_color="red")
             print(f"Load Error: {e}")
             return None

    def create_telemetry_page(self):
        ctk.CTkLabel(self.content_area, text="Telemetry Battle (Head-to-Head)", font=("Roboto", 24, "bold")).pack(pady=10)
        
        controls = ctk.CTkFrame(self.content_area)
        controls.pack(pady=10)
        
        self.d1_var = ctk.StringVar(value="VER")
        self.d2_var = ctk.StringVar(value="HAM")
        
        ctk.CTkEntry(controls, textvariable=self.d1_var, width=60).pack(side="left", padx=10)
        ctk.CTkLabel(controls, text="VS").pack(side="left", padx=10)
        ctk.CTkEntry(controls, textvariable=self.d2_var, width=60).pack(side="left", padx=10)
        
        ctk.CTkButton(controls, text="Analyze Fastest Lap", command=self.plot_telemetry, fg_color="#E10600").pack(side="left", padx=20)
        
        self.telemetry_frame = ctk.CTkFrame(self.content_area, fg_color="transparent")
        self.telemetry_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
    def plot_telemetry(self):
        session = self.load_session()
        if not session: return
        
        d1 = self.d1_var.get().upper()
        d2 = self.d2_var.get().upper()
        
        laps_d1 = session.laps.pick_drivers(d1).pick_fastest()
        laps_d2 = session.laps.pick_drivers(d2).pick_fastest()
        
        tel_d1 = laps_d1.get_car_data().add_distance()
        tel_d2 = laps_d2.get_car_data().add_distance()
        
        c1 = fastf1.plotting.get_driver_color(d1, session=session)
        c2 = fastf1.plotting.get_driver_color(d2, session=session)
        
        # Plot
        for widget in self.telemetry_frame.winfo_children(): widget.destroy()
        
        fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True, facecolor='#121212')
        
        # Speed
        ax[0].set_facecolor('#121212')
        ax[0].plot(tel_d1['Distance'], tel_d1['Speed'], color=c1, label=d1)
        ax[0].plot(tel_d2['Distance'], tel_d2['Speed'], color=c2, label=d2)
        ax[0].set_ylabel("Speed (km/h)", color="white")
        ax[0].tick_params(colors="white")
        ax[0].legend()
        ax[0].set_title(f"Speed Trace: {d1} vs {d2}", color="white")
        
        # Throttle/Brake
        ax[1].set_facecolor('#121212')
        ax[1].plot(tel_d1['Distance'], tel_d1['Throttle'], color=c1, linestyle="--", label=f"{d1} Throttle")
        ax[1].plot(tel_d2['Distance'], tel_d2['Throttle'], color=c2, linestyle="--", label=f"{d2} Throttle")
        
        # Overlay Brake? or separate? Let's just do Throttle for clarity or Brake as negative?
        # Requirement: "Throttle & Brake Inputs"
        # Let's add Brake on twin axis if possible, or just plot brake as well.
        
        ax[1].set_ylabel("Throttle %", color="white")
        ax[1].tick_params(colors="white")
        ax[1].set_xlabel("Distance (m)", color="white")
        
        canvas = FigureCanvasTkAgg(fig, master=self.telemetry_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas.draw()

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
                
                color = fastf1.plotting.COMPOUND_COLORS.get(compound, "white")
                
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
        ctk.CTkButton(btn_frame, text="Gear Shift Map", command=self.plot_gear_map, fg_color="#E10600").pack(side="left", padx=10)
        ctk.CTkButton(btn_frame, text="Speed Map", command=self.plot_speed_map, fg_color="#2CC985").pack(side="left", padx=10)
        
        self.track_frame = ctk.CTkFrame(self.content_area, fg_color="transparent")
        self.track_frame.pack(fill="both", expand=True, padx=10, pady=10)

    def plot_gear_map(self):
        session = self.load_session()
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
        session = self.load_session()
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
        for drv in drivers:
             laps = session.laps.pick_drivers(drv)
             try:
                 color = fastf1.plotting.get_driver_color(drv, session=session)
             except:
                 color = "white"
                 
             ax.plot(laps['LapNumber'], laps['Position'], color=color, label=drv)
             
        ax.set_ylabel("Position", color="white")
        ax.set_xlabel("Lap", color="white")
        ax.invert_yaxis()
        ax.tick_params(colors="white")
        ax.set_title("Position Changes", color="white")
        # ax.legend() # Too crowded
        
        canvas = FigureCanvasTkAgg(fig, master=self.pos_frame)
        canvas.get_tk_widget().pack(fill="both", expand=True)
        canvas.draw()

    def add_placeholder_plot(self, title):
        fig, ax = plt.subplots(figsize=(8, 5), facecolor='#121212')
        ax.set_facecolor('#121212')
        ax.text(0.5, 0.5, f"[ {title} Graph Here ]\n\nLoad Data to View", 
                color='white', ha='center', va='center', fontsize=14)
        ax.axis('off')
        
        canvas = FigureCanvasTkAgg(fig, master=self.content_area)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=20, pady=20)
        canvas.draw()

if __name__ == "__main__":
    app = AnalyticsDashboardApp()
    app.mainloop()