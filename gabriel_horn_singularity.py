import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button, CheckButtons
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp

def gabriel_horn(x, phi, scale=1.0):
    """
    Generate a Gabriel's Horn shape (Torricelli's Trumpet).
    
    Gabriel's Horn Paradox:
    - The surface is created by rotating y = 1/x around the x-axis
    - It has FINITE VOLUME but INFINITE SURFACE AREA
    - Volume = π (as x → ∞)
    - Surface Area = ∞ (diverges logarithmically)
    
    This paradox mirrors black hole singularities:
    - As x → 0 (singularity), the radius → ∞ (infinite curvature)
    - As x → ∞ (far from black hole), the radius → 0 (flat spacetime)
    
    Parameters:
    - x: x-coordinates (distance from center)
    - phi: angle in radians (0 to 2π)
    - scale: factor to scale the radius
    """
    # Prevent division by zero
    x_safe = np.maximum(x, 1e-10)
    r = scale / x_safe
    y = r * np.cos(phi)
    z = r * np.sin(phi)
    return x, y, z

def schwarzschild_metric(r, rs):
    """
    Calculate the Schwarzschild metric factor for gravitational redshift.
    
    This factor represents how time dilation and light redshift vary with distance
    from a black hole. The formula comes from Einstein's General Relativity:
    
    sqrt(1 - rs/r) where:
    - r: radial distance from black hole center
    - rs: Schwarzschild radius (event horizon)
    
    Returns:
    - 0 at or inside event horizon (infinite redshift)
    - Approaches 1 far from the black hole (no redshift)
    
    COLOR SIGNIFICANCE:
    - Red (factor → 0): Near event horizon, extreme gravitational redshift
    - Blue (factor → 1): Far from black hole, minimal redshift
    """
    # Avoid division by zero at the event horizon
    if np.isscalar(r):
        if r <= rs:
            return 0.0
        return np.sqrt(1.0 - rs/r)
    else:
        # Vectorized calculation for efficiency
        result = np.zeros_like(r)
        mask = r > rs
        result[mask] = np.sqrt(1.0 - rs/r[mask])
        return result

def photon_geodesic_schwarzschild(t, y, rs):
    """
    Calculate photon geodesics in Schwarzschild spacetime (simplified 2D).
    
    This simulates how light rays bend around a black hole.
    Uses simplified equations for visualization purposes.
    
    Parameters:
    - t: affine parameter (similar to time)
    - y: [r, phi, dr/dt, dphi/dt] - position and velocity in polar coordinates
    - rs: Schwarzschild radius
    
    Returns derivative vector for numerical integration.
    """
    r, phi, dr_dt, dphi_dt = y
    
    # Prevent singularity issues
    if r <= rs * 1.01:
        return [0, 0, 0, 0]
    
    # Simplified geodesic equations for photons in Schwarzschild metric
    # These are approximations for visualization
    f = 1 - rs / r
    
    # Second derivatives (acceleration due to gravity)
    d2r_dt2 = -rs / (2 * r**2) * dphi_dt**2 * r**2 + rs / (2 * r**2) * f * dr_dt**2
    d2phi_dt2 = -2 / r * dr_dt * dphi_dt
    
    return [dr_dt, dphi_dt, d2r_dt2, d2phi_dt2]

def calculate_light_ray(start_r, start_angle, impact_parameter, rs, max_t=50, z_offset=0.0, z_velocity=0.0):
    """
    Calculate a light ray trajectory around the black hole in 3D.
    
    Parameters:
    - start_r: starting radial distance
    - start_angle: starting angle (radians)
    - impact_parameter: determines how close the ray passes to black hole
    - rs: Schwarzschild radius
    - max_t: maximum integration time
    - z_offset: initial z-coordinate (for 3D trajectories)
    - z_velocity: initial z-velocity component
    
    Returns:
    - x, y, z coordinates of the light ray path
    """
    # Initial conditions: position and velocity
    r0 = start_r
    phi0 = start_angle
    dr0 = -0.5  # Radial velocity component (moving toward black hole)
    dphi0 = impact_parameter / r0**2  # Angular velocity
    
    y0 = [r0, phi0, dr0, dphi0]
    
    # Integrate the geodesic equations
    t_span = (0, max_t)
    t_eval = np.linspace(0, max_t, 200)
    
    try:
        sol = solve_ivp(
            lambda t, y: photon_geodesic_schwarzschild(t, y, rs),
            t_span, y0, t_eval=t_eval, method='RK45',
            events=lambda t, y: y[0] - rs * 1.01  # Stop at event horizon
        )
        
        r = sol.y[0]
        phi = sol.y[1]
        
        # Convert polar to Cartesian with 3D z-component
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        
        # Add 3D z-component that varies with trajectory
        # Simple model: z oscillates/decays as photon spirals
        t_normalized = sol.t / max_t
        z = z_offset * np.exp(-t_normalized * 0.5) * np.cos(phi * 2 + z_velocity)
        
        return x, y, z
    except:
        # If integration fails, return empty arrays
        return np.array([]), np.array([]), np.array([])

def black_hole_visualization():
    # Create the figure and 3D axis with adjusted spacing
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Initial parameters - optimized for visual clarity
    initial_schwarzschild_radius = 2.0  # Schwarzschild radius (event horizon)
    x_min, x_max = 0.1, 15.0  # Extended range for paradox visualization
    phi_resolution = 80  # Optimized resolution
    x_resolution = 120   # Optimized resolution
    
    # Create meshgrid for Gabriel's Horn
    x_coords = np.linspace(x_min, x_max, x_resolution)
    phi_coords = np.linspace(0, 2 * np.pi, phi_resolution)
    X, PHI = np.meshgrid(x_coords, phi_coords)
    
    # Event horizon sphere meshgrid
    theta_coords = np.linspace(0, np.pi, 50)
    phi_eh_coords = np.linspace(0, 2 * np.pi, 50)
    THETA, PHI_EH = np.meshgrid(theta_coords, phi_eh_coords)
    
    # Animation state
    animation_state = {
        'frame': 0,
        'running': False,
        'auto_rotate': False,
        'light_rays_data': [],  # Store pre-calculated light ray paths
        'animated_objects': []  # Store animated light ray segments
    }
    
    # Light ray visualization toggle
    show_light_rays = [True]  # Using list to make it mutable in nested function
    
    def update_event_horizon(schwarzschild_radius):
        """Generate the event horizon sphere coordinates."""
        x_eh = schwarzschild_radius * np.sin(THETA) * np.cos(PHI_EH)
        y_eh = schwarzschild_radius * np.sin(THETA) * np.sin(PHI_EH)
        z_eh = schwarzschild_radius * np.cos(THETA)
        return x_eh, y_eh, z_eh

    def update_photon_sphere(schwarzschild_radius):
        """Generate the photon sphere coordinates (1.5x Schwarzschild radius)."""
        ps_radius = 1.5 * schwarzschild_radius
        x_ps = ps_radius * np.sin(THETA) * np.cos(PHI_EH)
        y_ps = ps_radius * np.sin(THETA) * np.sin(PHI_EH)
        z_ps = ps_radius * np.cos(THETA)
        return x_ps, y_ps, z_ps
    
    def calculate_redshift(X, Y, Z, schwarzschild_radius):
        """
        Calculate gravitational redshift factor for each point.
        
        VECTORIZED for efficiency (avoids nested Python loops).
        
        The redshift factor determines the color:
        - 0.0 (Red): Maximum redshift, near event horizon
        - 1.0 (Blue): No redshift, far from black hole
        - Intermediate (Purple/Pink): Transition zone
        """
        # Vectorized calculation - much faster than nested loops
        r = np.sqrt(X**2 + Y**2 + Z**2)
        return schwarzschild_metric(r, schwarzschild_radius)
    
    def precalculate_light_rays(schwarzschild_radius):
        """Pre-calculate all light ray trajectories for animation."""
        animation_state['light_rays_data'] = []
        
        # Specific impact parameter multipliers to show different behaviors
        # Critical value is ~2.6 (3*sqrt(3)/2)
        # We want rays that get VERY close to the photon sphere (1.5 Rs)
        impact_multipliers = [
            1.8, 2.2, 2.4, 2.5,      # Captured (Red)
            2.58, 2.62, 2.7, 2.85,   # Critical/Strong Lensing (Orange)
            3.2, 4.0, 5.5, 7.0       # Weakly Deflected (Yellow)
        ]
        
        num_rays = len(impact_multipliers)
        start_radius = 15.0
        
        for i, mult in enumerate(impact_multipliers):
            angle = 2 * np.pi * i / num_rays
            impact_param = schwarzschild_radius * mult
            
            # Add 3D variation
            z_offset = (i % 3 - 1) * 2.0  # -2, 0, or +2
            z_velocity = i * 0.3
            
            x, y, z = calculate_light_ray(
                start_radius, angle, impact_param, 
                schwarzschild_radius, max_t=40,  # Increased time for spiraling
                z_offset=z_offset, z_velocity=z_velocity
            )
            
            if len(x) > 0:
                # Determine color based on trajectory
                min_r = np.min(np.sqrt(x**2 + y**2 + z**2))
                
                # Color logic based on user request:
                # Red: Captured
                # Orange: Critical scattering (near photon sphere at 1.5 Rs)
                # Yellow: Escaped
                
                if min_r < schwarzschild_radius * 1.05:
                    color = 'red'      # Captured
                    alpha = 0.9
                    linewidth = 2.5
                elif min_r < schwarzschild_radius * 1.8:
                    color = 'orange'   # Critical Scattering / Strong Lensing
                    alpha = 0.85
                    linewidth = 2.2
                else:
                    color = 'yellow'   # Escaped
                    alpha = 0.7
                    linewidth = 2.0
                
                animation_state['light_rays_data'].append({
                    'x': x, 'y': y, 'z': z,
                    'color': color, 'alpha': alpha, 'linewidth': linewidth
                })
    
    # Add color bar for redshift with detailed labels
    cbar_ax = fig.add_axes([0.92, 0.22, 0.015, 0.65])
    cbar = plt.colorbar(cm.ScalarMappable(cmap=cm.coolwarm), cax=cbar_ax)
    cbar.set_label('Gravitational Redshift Factor\n(0=Extreme, 1=None)', 
                   fontsize=8, labelpad=10)
    
    # Add legend text for light rays
    light_explanation = (
        "LIGHT RAY VISUALIZATION\n"
        "━━━━━━━━━━━━━━━━━━━━━━━\n"
        "🟡 YELLOW: Escaped\n"
        "   - Photons escape gravity\n\n"
        "🟠 ORANGE: Critical\n"
        "   - Strong lensing\n"
        "   - Near Photon Sphere\n\n"
        "🔴 RED: Captured\n"
        "   - Fall into Black Hole\n\n"
        "🔵 CYAN SPHERE:\n"
        "   - Photon Sphere (1.5 Rs)\n\n"
        "⭕ BRIGHT MARKERS:\n"
        "   - Real-time position"
    )
    
    # Add color explanation text
    color_explanation = (
        "SPACETIME COLOR:\n"
        "━━━━━━━━━━━━━━━━━━\n"
        "🔴 RED: Near horizon\n"
        "   - Extreme redshift\n"
        "   - Time dilation\n\n"
        "🟣 PURPLE: Transition\n"
        "   - Moderate effects\n\n"
        "🔵 BLUE: Far away\n"
        "   - Flat spacetime\n\n"
        "⚫ BLACK SPHERE:\n"
        "   Event horizon"
    )
    
    # Add text boxes
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.92)
    text_box1 = fig.text(0.02, 0.55, color_explanation, fontsize=7, 
             verticalalignment='center', bbox=props, family='monospace')
    
    props2 = dict(boxstyle='round', facecolor='lightcyan', alpha=0.92)
    text_box2 = fig.text(0.02, 0.20, light_explanation, fontsize=7,
             verticalalignment='center', bbox=props2, family='monospace')
    
    # Main plot updating function
    def update_plot(val):
        ax.clear()
        
        # Get slider values
        schwarzschild_radius = radius_slider.val
        scale_factor = scale_slider.val
        
        # Generate Gabriel's Horn geometry
        X_mesh, Y_mesh, Z_mesh = gabriel_horn(X, PHI, scale=scale_factor)
        
        # Calculate redshift (vectorized - efficient)
        redshift = calculate_redshift(X_mesh, Y_mesh, Z_mesh, schwarzschild_radius)
        
        # Plot Gabriel's Horn with gravitational redshift color mapping
        surf = ax.plot_surface(X_mesh, Y_mesh, Z_mesh, 
                             facecolors=cm.coolwarm(redshift), 
                             alpha=0.7, shade=True, linewidth=0, 
                             antialiased=True, edgecolor='none')
        
        # Plot event horizon sphere (Black)
        x_eh, y_eh, z_eh = update_event_horizon(schwarzschild_radius)
        horizon = ax.plot_surface(x_eh, y_eh, z_eh, color='black', 
                                alpha=0.95, shade=True, linewidth=0, zorder=10)

        # Plot Photon Sphere (Cyan, wireframe/transparent)
        x_ps, y_ps, z_ps = update_photon_sphere(schwarzschild_radius)
        photon_sphere = ax.plot_wireframe(x_ps, y_ps, z_ps, color='cyan', 
                                        alpha=0.15, linewidth=0.5, zorder=5)
        
        # Draw static light rays (full paths)
        if show_light_rays[0] and not animation_state['running']:
            for ray_data in animation_state['light_rays_data']:
                ax.plot(ray_data['x'], ray_data['y'], ray_data['z'],
                       color=ray_data['color'], alpha=ray_data['alpha']*0.5,
                       linewidth=ray_data['linewidth']*0.7, linestyle='-')
        
        # Add labels
        ax.text(0, 0, 0, "⚫\nSingularity", color='white', 
               fontsize=10, ha='center', weight='bold')
        ax.text(schwarzschild_radius, 0, schwarzschild_radius/2, 
               "Event Horizon", color='white', fontsize=8, ha='center')
        ax.text(schwarzschild_radius*1.5, 0, schwarzschild_radius*1.5, 
               "Photon Sphere", color='cyan', fontsize=8, ha='center')
        
        # Set plot limits with zoom control
        zoom_factor = zoom_slider.val
        ax.set_xlim(-x_max * zoom_factor, x_max * zoom_factor)
        ax.set_ylim(-x_max * zoom_factor, x_max * zoom_factor)
        ax.set_zlim(-x_max/2 * zoom_factor, x_max/2 * zoom_factor)
        
        # Labels
        ax.set_xlabel('X (Radial)', fontsize=9, labelpad=5)
        ax.set_ylabel('Y (Radial)', fontsize=9, labelpad=5)
        ax.set_zlabel('Z (Vertical)', fontsize=9, labelpad=5)
        
        # Title
        title_text = (
            'Black Hole Spacetime with Animated Light Ray Telemetry\n'
            'Gabriel\'s Horn Paradox: Finite Volume, Infinite Surface Area'
        )
        ax.set_title(title_text, fontsize=11, pad=12, weight='bold')
        
        # Update view angle (with auto-rotate if enabled)
        if animation_state['auto_rotate']:
            azim = (azimuth_slider.val + animation_state['frame']) % 360
        else:
            azim = azimuth_slider.val
        ax.view_init(elev=elevation_slider.val, azim=azim)
        
        # Visual enhancements
        ax.grid(True, alpha=0.2)
        ax.set_facecolor('#f5f5f5')
        
        # Update color bar
        cbar.update_normal(cm.ScalarMappable(cmap=cm.coolwarm))
        
        fig.canvas.draw_idle()
    
    # Animation function
    def animate(frame):
        if not animation_state['running']:
            return
        
        animation_state['frame'] = frame
        
        # Clear previous animated objects
        for obj in animation_state['animated_objects']:
            try:
                obj.remove()
            except:
                pass
        animation_state['animated_objects'].clear()
        
        # Update view if auto-rotate is enabled
        if animation_state['auto_rotate']:
            azim = (azimuth_slider.val + frame) % 360
            ax.view_init(elev=elevation_slider.val, azim=azim)
        
        # Animate light rays
        if show_light_rays[0]:
            for ray_data in animation_state['light_rays_data']:
                x, y, z = ray_data['x'], ray_data['y'], ray_data['z']
                total_points = len(x)
                
                if total_points > 0:
                    # Calculate how much of the path to show (progressive reveal)
                    reveal_length = 30  # Number of points in the "tail"
                    end_idx = (frame * 2) % (total_points + 50)  # Cycle through path
                    
                    if end_idx < total_points:
                        start_idx = max(0, end_idx - reveal_length)
                        
                        # Draw the animated segment
                        line = ax.plot(x[start_idx:end_idx], 
                                      y[start_idx:end_idx], 
                                      z[start_idx:end_idx],
                                      color=ray_data['color'], 
                                      alpha=ray_data['alpha'],
                                      linewidth=ray_data['linewidth'], 
                                      linestyle='-')[0]
                        animation_state['animated_objects'].append(line)
                        
                        # Draw a bright "photon" marker at the front
                        if end_idx > 0 and end_idx < total_points:
                            # Bright marker with white center and colored edge
                            marker = ax.scatter([x[end_idx-1]], [y[end_idx-1]], [z[end_idx-1]],
                                              c='white', s=120, 
                                              alpha=1.0, marker='o', edgecolors=ray_data['color'],
                                              linewidths=2.0, zorder=20)
                            animation_state['animated_objects'].append(marker)
        
        return animation_state['animated_objects']
    
    # Checkbox controls
    rax = plt.axes([0.02, 0.87, 0.13, 0.10], facecolor='lightgray')
    check = CheckButtons(rax, ['Show Light Rays', 'Auto-Rotate'], [True, False])
    
    def toggle_options(label):
        if label == 'Show Light Rays':
            show_light_rays[0] = not show_light_rays[0]
            update_plot(None)
        elif label == 'Auto-Rotate':
            animation_state['auto_rotate'] = not animation_state['auto_rotate']
    
    check.on_clicked(toggle_options)
    
    # Add sliders for interactive control
    slider_color = 'lightgoldenrodyellow'
    
    ax_radius = plt.axes([0.25, 0.18, 0.45, 0.02], facecolor=slider_color)
    radius_slider = Slider(
        ax=ax_radius,
        label='Event Horizon Radius',
        valmin=0.5,
        valmax=5.0,
        valinit=initial_schwarzschild_radius,
        valstep=0.1
    )
    
    ax_scale = plt.axes([0.25, 0.14, 0.45, 0.02], facecolor=slider_color)
    scale_slider = Slider(
        ax=ax_scale,
        label='Spacetime Curvature',
        valmin=0.1,
        valmax=2.0,
        valinit=1.0,
        valstep=0.05
    )
    
    ax_azimuth = plt.axes([0.25, 0.10, 0.45, 0.02], facecolor=slider_color)
    azimuth_slider = Slider(
        ax=ax_azimuth,
        label='Azimuth (Rotation)',
        valmin=0,
        valmax=360,
        valinit=45,
        valstep=5
    )
    
    ax_elevation = plt.axes([0.25, 0.06, 0.45, 0.02], facecolor=slider_color)
    elevation_slider = Slider(
        ax=ax_elevation,
        label='Elevation (Tilt)',
        valmin=-90,
        valmax=90,
        valinit=20,
        valstep=5
    )
    
    ax_zoom = plt.axes([0.25, 0.02, 0.45, 0.02], facecolor=slider_color)
    zoom_slider = Slider(
        ax=ax_zoom,
        label='Zoom',
        valmin=0.3,
        valmax=2.0,
        valinit=1.0,
        valstep=0.1
    )
    
    # Slider change handler
    def on_slider_change(val):
        # Recalculate light rays when schwarzschild radius changes
        schwarzschild_radius = radius_slider.val
        precalculate_light_rays(schwarzschild_radius)
        update_plot(val)
    
    # Connect sliders
    radius_slider.on_changed(on_slider_change)
    scale_slider.on_changed(update_plot)
    azimuth_slider.on_changed(update_plot)
    elevation_slider.on_changed(update_plot)
    zoom_slider.on_changed(update_plot)
    
    # Animation control buttons
    play_ax = plt.axes([0.72, 0.08, 0.06, 0.04])
    play_button = Button(play_ax, '▶️ Play', hovercolor='lightgreen')
    
    pause_ax = plt.axes([0.72, 0.03, 0.06, 0.04])
    pause_button = Button(pause_ax, '⏸️ Pause', hovercolor='lightcoral')
    
    reset_ax = plt.axes([0.79, 0.055, 0.06, 0.04])
    reset_button = Button(reset_ax, '🔄 Reset', hovercolor='lightblue')
    
    def play_animation(event):
        animation_state['running'] = True
        animation_state['frame'] = 0
    
    def pause_animation(event):
        animation_state['running'] = False
        update_plot(None)  # Redraw static view
    
    def reset(event):
        """Reset all sliders to default values."""
        animation_state['running'] = False
        animation_state['frame'] = 0
        radius_slider.reset()
        scale_slider.reset()
        azimuth_slider.reset()
        elevation_slider.reset()
        zoom_slider.reset()
    
    play_button.on_clicked(play_animation)
    pause_button.on_clicked(pause_animation)
    reset_button.on_clicked(reset)
    
    # Pre-calculate light rays for initial view
    precalculate_light_rays(initial_schwarzschild_radius)
    
    # Initial plot
    update_plot(None)
    
    # Create animation object
    anim = FuncAnimation(fig, animate, frames=300, interval=50, 
                        blit=False, repeat=True, cache_frame_data=False)
    
    # Adjust layout
    plt.subplots_adjust(left=0.16, right=0.90, bottom=0.22, top=0.94)
    plt.show()

if __name__ == "__main__":
    # Force UTF-8 encoding for Windows console
    import sys
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

    print("=" * 80)
    print("ANIMATED BLACK HOLE VISUALIZATION WITH 3D PHOTON TELEMETRY")
    print("=" * 80)
    print("\n*** ANIMATION FEATURES:")
    print("-" * 80)
    print("PLAY button   : Start animating light rays traveling through spacetime")
    print("PAUSE button  : Pause animation and show complete light ray paths")
    print("RESET button  : Reset all parameters to default values")
    print()
    print("Show Light Rays: Toggle light ray visualization on/off")
    print("Auto-Rotate    : Enable automatic camera rotation during animation")
    print("-" * 80)
    print("\n*** LIGHT RAY ANIMATION (Now in 3D!):")
    print("-" * 80)
    print("Watch photons (light particles) travel along their curved paths!")
    print()
    print("YELLOW rays : Escaping photons")
    print("ORANGE rays : Critical scattering (near Photon Sphere)")
    print("RED rays    : Captured photons")
    print("CYAN sphere : Photon Sphere (1.5x Schwarzschild Radius)")
    print("=" * 80)
    print("\n*** COLOR SCHEME:")
    print("RED spacetime   : Extreme gravitational redshift (near horizon)")
    print("PURPLE spacetime: Medium gravitational effects")
    print("BLUE spacetime  : Weak gravitational effects (far away)")
    print("BLACK sphere   : Event horizon (point of no return)")
    print("-" * 80)
    print("\n*** THE PARADOX:")
    print("Gabriel's Horn: FINITE VOLUME but INFINITE SURFACE AREA")
    print("Black Hole: Finite mass but INFINITE density at singularity")
    print("=" * 80)
    print("\n*** READY TO EXPLORE!")
    print("Click 'Play' to watch light rays travel through curved spacetime!")
    print("Enable 'Auto-Rotate' for a dynamic 360-degree view!\n")
    
    black_hole_visualization()