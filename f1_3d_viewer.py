import sys
import json
import logging
import threading
import signal
import time
import numpy as np
from vispy import scene, app, io
from vispy.scene.cameras.turntable import TurntableCamera

# Setup simple logging to stdout instead of a file
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
TRACK_COLOR = "#333333"
CURTAIN_COLOR = "#1C1C1C"
BG_COLOR = "#121212"

def hex_to_rgba(hex_code, alpha=1.0):
    hex_code = hex_code.lstrip('#')
    return tuple(int(hex_code[i:i+2], 16)/255.0 for i in (0, 2, 4)) + (alpha,)

# -----------------------------------------------------------------------------
# GLOBAL STATE
# -----------------------------------------------------------------------------
class AppState:
    def __init__(self):
        self.driver_dots = {}
        self.driver_labels = {}
        self.lock = threading.Lock()
        self.canvas = None
        self.view = None
        self.track_line = None
        self.curtain_line = None
        self.hud_text = None
        self.controls_text = None
        self.is_running = True
        self.target_frame_data = {} # Latest data from main app
        self.current_frame_data = {} # LERP interpolated data
        
        # UI State Variables
        self.speed_str = "1x"
        self.lap_time_str = "00:00.000"
        self.is_paused = False
        self.finish_status = ""

state = AppState()

# -----------------------------------------------------------------------------
# IPC THREAD (Reads from stdin)
# -----------------------------------------------------------------------------
def read_input_stream():
    """ Runs in a background thread to read commands from the main Tkinter app. """
    logging.info("Input stream reader started.")
    try:
        for line in sys.stdin:
            if not state.is_running:
                break
                
            line = line.strip()
            if not line: continue
            
            try:
                data = json.loads(line)
                
                cmd = data.get('cmd')
                
                if cmd == 'init_track':
                    # Received track geometry
                    x = np.array(data['x'], dtype=np.float32)
                    y = np.array(data['y'], dtype=np.float32)
                    z = np.array(data['z'], dtype=np.float32)
                    bounds = data.get('bounds')
                    
                    logging.info(f"Received init_track with {len(x)} points.")
                    # Update canvas safely via Qt/VisPy event system or just rely on the next frame loop picking it up.
                    # For simplicity, we flag it so the main loop builds it. 
                    # Actually, we can just build it here if we are careful, but safer in main thread.
                    with state.lock:
                        state.track_data = (x, y, z)
                        state.track_bounds = bounds
                        state.needs_track_build = True

                elif cmd == 'update_frame':
                    # Received driver positions and HUD info
                    drivers = data.get('drivers', {})
                    
                    with state.lock:
                        # Store as Target Data for LERP
                        state.target_frame_data = drivers
                        
                        # Update HUD values
                        if 'ui' in data:
                             state.speed_str = data['ui'].get('speed', state.speed_str)
                             state.lap_time_str = data['ui'].get('lap_time', state.lap_time_str)
                             state.is_paused = data['ui'].get('is_paused', state.is_paused)
                             state.finish_status = data['ui'].get('finish_status', state.finish_status)
                             
                             # If paused, snap current to target immediately to prevent drifting after pause
                             if getattr(state, 'was_paused', False) and not state.is_paused:
                                  state.current_frame_data = {k: v.copy() for k, v in state.target_frame_data.items()}
                             state.was_paused = state.is_paused

                elif cmd == 'exit':
                    logging.info("Received exit command.")
                    state.is_running = False
                    app.quit()
                    break
                    
            except json.JSONDecodeError as e:
                logging.error(f"JSON Parse Error: {e} - Line: {line[:50]}...")
            except Exception as e:
                logging.error(f"Error processing command: {e}")
                
    except Exception as e:
         logging.error(f"Fatal Input Stream Error: {e}")
    finally:
         logging.info("Input stream reader thread exiting.")

# -----------------------------------------------------------------------------
# CUSTOM INERTIA CAMERA
# -----------------------------------------------------------------------------
class InertiaTurntableCamera(TurntableCamera):
    """ Custom camera that overrides mouse handling to add inertia perfectly matching Matplotlib """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Disable default interactive to use our own event handlers
        self.interactive = False
        
        self.drag_start = None
        self.azimuth_vel = 0.5  # Start with gentle auto-rotation
        self.elevation_vel = 0.0
        self.friction = 0.85 # Reduced friction so it stops faster (less floaty)
        
    def viewbox_mouse_event(self, event):
        # We only care about left click (button=1) and scroll (type='mouse_wheel')
        # We must intercept mouse_move events that happen while dragging
        
        if event.type == 'mouse_wheel':
            # Let the default handle zoom
            super().viewbox_mouse_event(event)
            return

        if event.button == 1:
            if event.type == 'mouse_press':
                self.drag_start = (event.pos[0], event.pos[1])
                self.azimuth_vel = 0.0
                self.elevation_vel = 0.0
                event.handled = True
                
            elif event.type == 'mouse_move' and self.drag_start is not None:
                dx = event.pos[0] - self.drag_start[0]
                dy = event.pos[1] - self.drag_start[1]
                
                # Inverted tilt direction and smoothed velocity for better inertia release
                new_azi = dx * -0.2 
                new_ele = dy * 0.2
                
                self.azimuth_vel = self.azimuth_vel * 0.5 + new_azi * 0.5
                self.elevation_vel = self.elevation_vel * 0.5 + new_ele * 0.5
                
                self.azimuth += self.azimuth_vel
                self.elevation = np.clip(self.elevation + self.elevation_vel, 0, 90)
                
                self.drag_start = (event.pos[0], event.pos[1])
                self.view_changed()
                event.handled = True
                
            elif event.type == 'mouse_release':
                self.drag_start = None
                event.handled = True
        
        # Explicitly skip passing generic mouse_move events up to base
        # This prevents the "Scale" and "is_dragging" attribute errors
        # arising from VisPy's default handler trying to process our intercepted states.

    def apply_inertia(self):
        if self.drag_start is None: # Not dragging
            if abs(self.azimuth_vel) > 0.01 or abs(self.elevation_vel) > 0.01:
                self.azimuth += self.azimuth_vel
                self.elevation = np.clip(self.elevation + self.elevation_vel, 0, 90)
                
                self.azimuth_vel *= self.friction
                self.elevation_vel *= self.friction
                
                self.view_changed()


# -----------------------------------------------------------------------------
# VISPY MAIN LOOP UPDATE (60 FPS)
# -----------------------------------------------------------------------------
def on_timer(event):
    if not state.is_running:
        return
        
    # Apply Camera Physics
    if state.view and isinstance(state.view.camera, InertiaTurntableCamera):
        state.view.camera.apply_inertia()

    with state.lock:
        # 1. Build Track if needed
        if getattr(state, 'needs_track_build', False) and state.view:
            x, y, z = state.track_data
            
            track_points = np.column_stack((x, y, z))
            state.track_line = scene.visuals.Line(pos=track_points, 
                                            color=hex_to_rgba(TRACK_COLOR, 0.5), 
                                            width=4, 
                                            parent=state.view.scene)
            
            z_floor = np.min(z) - 20
            floor_points = np.column_stack((x, y, np.full_like(z, z_floor)))
            
            curtain_points = np.zeros((2 * len(x), 3), dtype=np.float32)
            curtain_points[0::2] = track_points
            curtain_points[1::2] = floor_points
            
            state.curtain_line = scene.visuals.Line(pos=curtain_points, 
                                              connect='segments', 
                                              color=hex_to_rgba(CURTAIN_COLOR, 0.4), 
                                              width=1, 
                                              parent=state.view.scene)
            
            # Setup Camera Bounds
            if hasattr(state, 'track_bounds') and state.track_bounds:
                min_x, max_x, min_y, max_y, min_z, max_z = state.track_bounds
                center = ((min_x+max_x)/2, (min_y+max_y)/2, (min_z+max_z)/2)
                state.view.camera.center = center
                # Approximate distance based on track size
                state.view.camera.distance = max(max_x-min_x, max_y-min_y) * 1.2
                state.view.camera.elevation = 45 # Default angle
            
            # Setup HUD
            state.hud_text = scene.visuals.Text("", color='white', font_size=16, bold=True, 
                                                anchor_x='left', anchor_y='top', parent=state.canvas.scene)
            # FIX: Moved down from [10,20] to ensure it clears the window borders
            state.hud_text.pos = [20, 40] 
            
            state.controls_text = scene.visuals.Text("CONTROLS:  [SPACE] Pause/Resume  |  [UP/DOWN] Replay Speed", 
                                                     color='#AAAAAA', font_size=10, 
                                                     anchor_x='center', anchor_y='bottom', parent=state.canvas.scene)
            # Position at bottom center
            state.controls_text.pos = [state.canvas.size[0] / 2, state.canvas.size[1] - 20]

            state.needs_track_build = False
        
        # Handle Resize for HUD
        if getattr(state, 'controls_text', None):
             state.controls_text.pos = [state.canvas.size[0] / 2, state.canvas.size[1] - 20]
        
        # Update HUD Text
        if getattr(state, 'hud_text', None):
             status = "PAUSED" if state.is_paused else ("FINISHED" if state.finish_status else state.speed_str)
             color = "#E10600" if state.is_paused else ("#AAAAAA" if state.finish_status else "#2CC985")
             
             display_str = f"Time: {state.lap_time_str}  |  {status}"
             if state.finish_status: display_str += f" - {state.finish_status}"
             
             state.hud_text.text = display_str
             state.hud_text.color = color

        
        # 2. Update Drivers (LERP Interpolation)
        target_data = getattr(state, 'target_frame_data', None)
        if target_data and state.view:
            # LERP Factor (Higher = faster snap. 0.2 at 60fps is smooth but responsive)
            # We bypass LERP if paused so dots don't drift after hitting pause
            lerp_t = 1.0 if state.is_paused else 0.2 
            
            for driver, info in target_data.items():
                tx, ty, tz = info['x'], info['y'], info['z']
                color = info['color']
                visible = info['visible']
                
                # Initialize current state if new driver
                if driver not in state.current_frame_data:
                     state.current_frame_data[driver] = {'x': tx, 'y': ty, 'z': tz}
                     
                cx = state.current_frame_data[driver]['x']
                cy = state.current_frame_data[driver]['y']
                cz = state.current_frame_data[driver]['z']
                
                # Handle large jumps (teleportation across map) -> snap instead of lerp
                dx = tx - cx
                dy = ty - cy
                if abs(dx) > 1000 or abs(dy) > 1000:
                     cx, cy, cz = tx, ty, tz
                else:
                     # Calculate LERP
                     cx = cx + (tx - cx) * lerp_t
                     cy = cy + (ty - cy) * lerp_t
                     cz = cz + (tz - cz) * lerp_t
                
                # Save interpolated back to state
                state.current_frame_data[driver]['x'] = cx
                state.current_frame_data[driver]['y'] = cy
                state.current_frame_data[driver]['z'] = cz
                
                
                if driver not in state.driver_dots:
                    # Create new marker and text
                    dot = scene.visuals.Markers(parent=state.view.scene)
                    text = scene.visuals.Text(text=driver, color=color, bold=True, font_size=10, parent=state.view.scene)
                    state.driver_dots[driver] = dot
                    state.driver_labels[driver] = text
                
                dot = state.driver_dots[driver]
                text = state.driver_labels[driver]
                
                if visible and not np.isnan(cx) and not np.isnan(cy):
                     # Elevate Z slightly so it doesn't clip into the track line
                     dot.set_data(pos=np.array([[cx, cy, cz + 1.5]]), face_color=hex_to_rgba(color), edge_color=(0,0,0,1), size=8)
                     # Position text slightly above
                     text.pos = [cx, cy, cz + 6.5]
                     
                     # Simple logic for DNF
                     if info.get('is_dnf', False):
                          dot.set_data(pos=np.array([[cx, cy, cz]]), face_color=hex_to_rgba("#666666"), edge_color=(0,0,0,1), size=8)
                          text.color = "#666666"
                          text.text = f"{driver} (DNF)"
                     else:
                          text.color = color
                          text.text = driver 
                else:
                     dot.set_data(pos=np.empty((0,3)))
                     text.text = ""

    state.canvas.update()

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
def main():
    try:
        # Check if PyQt5/VisPy is actually working before committing to the full loop
        # This will throw an error immediately if OpenGL is missing, allowing the parent to catch it.
        state.canvas = scene.SceneCanvas(keys='interactive', show=True, bgcolor=hex_to_rgba(BG_COLOR), title="F1 Race Replayer - 3D Engine")
    except Exception as e:
        # Print critical error to stderr so the parent process can read it immediately
        print(f"CRITICAL_VISPY_INIT_ERROR: {e}", file=sys.stderr)
        sys.exit(1)
        
    state.view = state.canvas.central_widget.add_view()
    # Use the Custom Inertia Camera
    state.view.camera = InertiaTurntableCamera()
    
    # Handle Keyboard Controls (IPC back to main app)
    @state.canvas.events.key_press.connect
    def on_key_press(event):
        key = event.key.name
        cmd = None
        if key == 'Space':
            cmd = {"cmd": "input", "action": "toggle_pause"}
        elif key == 'Up':
            cmd = {"cmd": "input", "action": "speed_up"}
        elif key == 'Down':
            cmd = {"cmd": "input", "action": "speed_down"}
            
        if cmd:
            # Send back to Tkinter via stdout. 
            # Tkinter needs to listen to this subprocess's stdout.
            print(json.dumps(cmd))
            sys.stdout.flush()
            
    # Resize handles UI
    @state.canvas.events.resize.connect
    def on_resize(event):
        if getattr(state, 'controls_text', None):
             state.controls_text.pos = [event.size[0] / 2, event.size[1] - 20]
    
    # Handle window close event
    @state.canvas.events.close.connect
    def on_close(event):
        logging.info("Window closed. Exiting.")
        state.is_running = False
        app.quit()
        sys.exit(0)

    # Start IPC Thread
    ipc_thread = threading.Thread(target=read_input_stream, daemon=True)
    ipc_thread.start()

    # Start 60FPS loop
    timer = app.Timer(interval=1/60.0, connect=on_timer, start=True)
    
    # Signal readiness to parent
    print("VISPY_READY")
    sys.stdout.flush()

    logging.info("Starting VisPy App Loop")
    if sys.flags.interactive != 1:
        app.run()

if __name__ == '__main__':
    main()
