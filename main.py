import customtkinter as ctk
import subprocess
import sys
import os

# Configure appearance
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

class MainApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("F1 Visualizer Hub")
        self.geometry("800x600")
        
        # Center the window
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        
        self.main_frame = ctk.CTkFrame(self, corner_radius=15)
        self.main_frame.grid(row=0, column=0, padx=40, pady=40, sticky="nsew")
        
        self.main_frame.grid_columnconfigure(0, weight=1)
        
        # Header
        ctk.CTkLabel(self.main_frame, text="F1 VISUALIZER", font=("Roboto", 40, "bold"), text_color="#E10600").pack(pady=(60, 20))
        ctk.CTkLabel(self.main_frame, text="Select Mode", font=("Roboto", 20), text_color="gray").pack(pady=(0, 40))
        
        # Buttons
        self.btn_replayer = ctk.CTkButton(self.main_frame, text="RACE REPLAYER", 
                                          command=self.launch_replayer,
                                          width=300, height=80,
                                          font=("Roboto", 20, "bold"),
                                          fg_color="#E10600", hover_color="#B00500")
        self.btn_replayer.pack(pady=20)
        
        self.btn_dashboard = ctk.CTkButton(self.main_frame, text="ANALYTICS DASHBOARD", 
                                           command=self.launch_dashboard,
                                           width=300, height=80,
                                           font=("Roboto", 20, "bold"),
                                           fg_color="#1f538d", hover_color="#14375e")
        self.btn_dashboard.pack(pady=20)
        
        # Footer
        ctk.CTkLabel(self.main_frame, text="Powered by FastF1", font=("Mono", 10), text_color="gray50").pack(side="bottom", pady=20)

    def launch_replayer(self):
        # Launch race_replayer.py as a separate process
        script_path = os.path.join(os.path.dirname(__file__), "race_replayer.py")
        print(f"Launching: {script_path}")
        subprocess.Popen([sys.executable, script_path])

    def launch_dashboard(self):
        # Launch analytics_dashboard.py as a separate process
        script_path = os.path.join(os.path.dirname(__file__), "analytics_dashboard.py")
        print(f"Launching: {script_path}")
        subprocess.Popen([sys.executable, script_path])

if __name__ == "__main__":
    app = MainApp()
    app.mainloop()
