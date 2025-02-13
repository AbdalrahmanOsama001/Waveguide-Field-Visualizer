import sys
import os
import gc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, font
import webbrowser

def scale_vectors(u, v, desired_length):
    """Normalize vector fields to fixed length while preserving direction.
    
    Args:
        u (ndarray): X-component of vector field
        v (ndarray): Y-component of vector field
        desired_length (float): Target length for all vectors
        
    Returns:
        tuple: (u_scaled, v_scaled) - Normalized vector components
        
    Note:
        Handles zero-length vectors safely by returning original components
    """
    magnitudes = np.sqrt(u**2 + v**2)
    scale_factors = np.ones_like(magnitudes) * desired_length
    valid_mask = magnitudes > 1e-10
    scale_factors[valid_mask] = desired_length / magnitudes[valid_mask]
    return u * scale_factors, v * scale_factors

class WaveguideApp:
    """Main application class for waveguide electromagnetic field visualization.
    
    Attributes:
        root (tk.Tk): Main Tkinter root window
        mode_var (tk.StringVar): TE/TM mode selection variable
        m_var (tk.StringVar): Mode index m input variable
        n_var (tk.StringVar): Mode index n input variable
        z_var (tk.StringVar): Z-position input variable
        y_cut_var (tk.StringVar): Y-cut position input variable
        x_cut_var (tk.StringVar): X-cut position input variable
        plot_frame (ttk.Frame): Frame container for matplotlib visualizations
    """

    def __init__(self, root):
        """Initialize the Waveguide Visualizer application GUI.
        
        Args:
            root (tk.Tk): Root window for the Tkinter application
        """
        self.root = root
        self.root.title("Waveguide Field Visualizer")
        
        # Input Frame
        input_frame = ttk.Frame(root, padding="10")
        input_frame.grid(row=0, column=0, sticky="nw")
        
        self.mode_var = tk.StringVar()
        self.m_var = tk.StringVar()
        self.n_var = tk.StringVar()
        self.z_var = tk.StringVar(value="0.0")
        self.y_cut_var = tk.StringVar(value="0.5")
        self.x_cut_var = tk.StringVar(value="0.5")
        
        ttk.Label(input_frame, text="Mode (1=TE, 2=TM):").grid(row=0, column=0, sticky="w")
        ttk.Entry(input_frame, textvariable=self.mode_var, width=5).grid(row=0, column=1)
        
        ttk.Label(input_frame, text="m:").grid(row=1, column=0, sticky="w")
        ttk.Entry(input_frame, textvariable=self.m_var, width=5).grid(row=1, column=1)
        
        ttk.Label(input_frame, text="n:").grid(row=2, column=0, sticky="w")
        ttk.Entry(input_frame, textvariable=self.n_var, width=5).grid(row=2, column=1)
        
        ttk.Label(input_frame, text="Z (λ):").grid(row=3, column=0, sticky="w")
        ttk.Entry(input_frame, textvariable=self.z_var, width=5).grid(row=3, column=1)
        
        ttk.Label(input_frame, text="Y-cut (0-1):").grid(row=4, column=0, sticky="w")
        ttk.Entry(input_frame, textvariable=self.y_cut_var, width=5).grid(row=4, column=1)
        
        ttk.Label(input_frame, text="X-cut (0-1):").grid(row=5, column=0, sticky="w")
        ttk.Entry(input_frame, textvariable=self.x_cut_var, width=5).grid(row=5, column=1)
        
        ttk.Button(input_frame, text="Generate Plot", command=self.generate_plot).grid(row=6, column=0, columnspan=2, pady=10)
        
        # Plot Frame
        self.plot_frame = ttk.Frame(root)
        self.plot_frame.grid(row=0, column=1, sticky="nsew")
        self.plot_frame.grid_rowconfigure(0, weight=1)
        self.plot_frame.grid_columnconfigure(0, weight=1)
        
        # Credits Frame
        self.credits_frame = ttk.LabelFrame(root, text="Credits", padding="10")
        self.credits_frame.grid(row=1, column=0, columnspan=2, sticky="ew")
        
        def open_url(url):
            webbrowser.open(url)
        
        credits_text1 = "Created by Abdalrahman Osama | CCE-E'24"
        github_link = "GitHub Repository"
        
        credits_label = tk.Label(self.credits_frame, text=credits_text1, fg="black", cursor="arrow", font="Helvetica 10 bold")
        credits_label.pack(side=tk.TOP)
        
        github_label = tk.Label(self.credits_frame, text=github_link, fg="blue", cursor="hand2", font="Helvetica 10 bold")
        github_label.pack(side=tk.TOP)
        github_label.bind("<Button-1>", lambda e: open_url("https://github.com/AbdalrahmanOsama001/Waveguide-Field-Visualizer"))
        
        # Add underline to hyperlink
        link_font = font.Font(github_label, github_label.cget("font"))
        link_font.configure(underline=True)
        github_label.configure(font=link_font)
        
        # Configure grid weights
        root.grid_columnconfigure(1, weight=1)
        root.grid_rowconfigure(0, weight=1)
        root.grid_rowconfigure(1, weight=0)
        
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.active_figures = []  # Track all figures to close on exit
        
    def on_close(self):
        """Force-clean all resources"""
        # Close all matplotlib figures
        plt.close('all')
        # Destroy all plot elements
        if hasattr(self, 'canvas'):
            self.canvas.get_tk_widget().destroy()
            del self.canvas
        # Force garbage collection
        gc.collect()
        # Nuclear option for process termination
        self.root.destroy()
        os._exit(0)  # Harshest termination

    def generate_plot(self):
        """Validate inputs and generate updated field visualizations.
        
        Performs input validation, clears previous plots, and triggers
        the creation of new field visualizations. Displays error messages
        for invalid inputs using tkinter messagebox.
        """
        try:
            mode = int(self.mode_var.get())
            m = int(self.m_var.get())
            n = int(self.n_var.get())
            z = float(self.z_var.get())
            y_cut = float(self.y_cut_var.get())
            x_cut = float(self.x_cut_var.get())
        except ValueError:
            tk.messagebox.showerror("Error", "Invalid input values")
            return

        # Clear previous plots
        for widget in self.plot_frame.winfo_children():
            widget.destroy()
        # Explicitly close previous figures
        for fig in self.active_figures:
            plt.close(fig)
        self.active_figures = []
        

        # Create new figure and canvas
        fig = plt.figure(figsize=(10, 8))
        self.active_figures.append(fig)
        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # Generate plots
        self.create_plots(fig, mode, m, n, z, y_cut, x_cut)
        canvas.draw()
        
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def create_plots(self, fig, mode, m, n, z, y_cut, x_cut):
        """Create the 2x2 plot grid with field visualizations.
        
        Args:
            fig (matplotlib.figure.Figure): Figure object for plotting
            mode (int): Propagation mode (1=TE, 2=TM)
            m (int): Mode index for x-direction
            n (int): Mode index for y-direction
            z (float): Z-position in wavelengths (λ)
            y_cut (float): Y-cut position (0-1 relative to height)
            x_cut (float): X-cut position (0-1 relative to width)
        """
        axs = fig.subplots(2, 2)
        MODE = "TE" if mode == 1 else "TM"
        fig.suptitle(f"Rectangular Waveguide Field Plots - Mode: {MODE}{m}{n}", fontsize=14)
        
        self.waveguide_transversal_section(mode, m, n, z, axs[0, 0])
        self.waveguide_side_section(mode, m, n, z, x_cut, axs[0, 1])
        self.waveguide_top_section(mode, m, n, z, y_cut, axs[1, 0])
        
        # Legend area now only contains the description text (credits removed)
        ax_legend = axs[1, 1]
        ax_legend.axis('off')
        description_text = """Color Description:
E-field: Blue
H-field: Red

Marker Description:
Arrows: In-plane components
'o': Out-of-plane negative
'x': In-plane positive"""
        ax_legend.text(0.5, 0.5, description_text, ha='center', va='center')
        
        plt.tight_layout(rect=[0, 0.05, 1, 0.90])

    def waveguide_transversal_section(self, mode, m, n, z, ax):
        """Plot transversal (XY-plane) field components at specified Z-position.
        
        Args:
            mode (int): 1=TE mode, 2=TM mode
            m (int): X-direction mode index
            n (int): Y-direction mode index
            z (float): Z-position in wavelengths (λ)
            ax (matplotlib.axes.Axes): Axes object for plotting
            
        Calculates:
            - E-field components (Ex, Ey) using sinusoidal mode patterns
            - H-field components (Hx, Hy, Hz) from Maxwell's equations
            - Vector fields normalized for consistent arrow lengths
        """
        B = 2 * np.pi
        a, b = 1, 1
        regions_x = 2 * m + 1
        regions_y = 2 * n + 1
        x_vals = np.linspace(0, a, regions_x)
        y_vals = np.linspace(0, b, regions_y)
        if n == 0:  # if n==0, center y at b/2
            y_vals = np.array([b/2])
        X, Y = np.meshgrid(x_vals, y_vals, indexing='ij')
        Ex, Ey, Hx, Hy, Hz, Ez = (np.zeros(X.shape) for _ in range(6))

        cosz = np.cos(B * z)
        sinz = np.sin(B * z)
        sinx = np.sin(m * np.pi * X / a)
        cosx = np.cos(m * np.pi * X / a)
        siny = np.sin(n * np.pi * Y / b)
        cosy = np.cos(n * np.pi * Y / b)

        if mode == 1:  # TE Mode
            Ex = cosx * siny * sinz
            Hy = Ex
            Ey = -sinx * cosy * sinz
            Hx = -Ey
            Hz = cosx * cosy * cosz
        elif mode == 2:  # TM Mode
            Ex = -cosx * siny * sinz
            Hy = Ex
            Ey = -sinx * cosy * sinz
            Hx = -Ey
            Ez = sinx * siny * cosz
        else:
            print("Invalid mode")
            return

        ax.set_xticks(x_vals)
        ax.set_yticks(y_vals)
        ax.grid(color='k', linestyle='--', linewidth=0.5)
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
        ax.set_title(f"Transversal Section (XY) at Z = {z:.2f}λ")
        
        grid_spacing_x = a / (regions_x - 1)
        grid_spacing_y = b if regions_y == 1 else b / (regions_y - 1)  # avoid division by zero when n==0
        arrow_length = 0.4 * min(grid_spacing_x, grid_spacing_y)
        
        Ex_scaled, Ey_scaled = scale_vectors(Ex, Ey, arrow_length)
        Hx_scaled, Hy_scaled = scale_vectors(Hx, Hy, arrow_length)
        
        ax.quiver(X, Y, Ex_scaled, Ey_scaled, 
                color='b', scale=1, scale_units='xy',
                width=0.006, headwidth=5)
        ax.quiver(X, Y, Hx_scaled, Hy_scaled,
                color='r', scale=1, scale_units='xy',
                width=0.006, headwidth=5)

    def waveguide_top_section(self, mode, m, n, z, y_cut, ax):
        """Plot top (XZ-plane) field components at fixed Y-cut position.
        
        Args:
            mode (int): 1=TE mode, 2=TM mode
            m (int): X-direction mode index
            n (int): Y-direction mode index
            z (float): Central Z-position in wavelengths (λ)
            y_cut (float): Y-cut position (0-1 relative to height)
            ax (matplotlib.axes.Axes): Axes object for plotting
            
        Features:
            - Shows longitudinal (Z-axis) field variations
            - Displays out-of-plane components using markers:
              * 'o' for negative Y-direction fields
              * 'x' for positive Y-direction fields
            - Inverts Y-axis to match waveguide propagation convention
        """
        B = 2 * np.pi
        a, b = 1, 1
        z_vals = np.array([z - 0.5, z - 0.25, z, z + 0.25, z + 0.5])
        Y_fixed = np.clip(y_cut * b, 0, b)
        if n == 0:  # if n==0, center y at b/2
            Y_fixed = b/2
        regions_x = 2 * m + 1
        x_vals = np.linspace(0, a, regions_x)
        X, Z = np.meshgrid(x_vals, z_vals, indexing='ij')

        sinx = np.sin(m * np.pi * X / a)
        cosx = np.cos(m * np.pi * X / a)
        siny_val = np.sin(n * np.pi * Y_fixed / b)
        cosy_val = np.cos(n * np.pi * Y_fixed / b)
        sinz = np.sin(B * Z)
        cosz = np.cos(B * Z)

        Ex_top, Ez_top = np.zeros_like(X), np.zeros_like(X)
        Hx_top, Hz_top = np.zeros_like(X), np.zeros_like(X)
        Ey_top, Hy_top = np.zeros_like(X), np.zeros_like(X)

        if mode == 1:  # TE Mode
            Hz_top = -cosx * cosy_val * cosz
            Ex_top = cosx * siny_val * sinz
            Hy_top = Ex_top
            Ey_top = -sinx * cosy_val * sinz
            Hx_top = -Ey_top

        elif mode == 2:  # TM Mode
            Ez_top = -sinx * siny_val * cosz
            Ex_top = -cosx * siny_val * sinz
            Ey_top = -sinx * cosy_val * sinz
            Hx_top = (n * np.pi / b) * sinx * cosy_val * sinz
            Hy_top = -(m * np.pi / a) * cosx * siny_val * sinz

        threshold = 1e-10
        Ex_top[np.abs(Ex_top) < threshold] = 0
        Ey_top[np.abs(Ey_top) < threshold] = 0
        Hx_top[np.abs(Hx_top) < threshold] = 0
        Hy_top[np.abs(Hy_top) < threshold] = 0

        ax.invert_yaxis()
        ax.set_xticks(x_vals)
        ax.set_yticks(z_vals)
        ax.grid(True, linestyle='--', linewidth=0.5)
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Z-axis (downwards)")
        ax.set_title(f"Top Section (XZ) at Y = {Y_fixed:.2f}b")

        grid_spacing_x = a / (regions_x - 1)
        arrow_length = 0.4 * grid_spacing_x
        
        Ex_scaled, Ez_scaled = scale_vectors(Ex_top, Ez_top, arrow_length)
        Hx_scaled, Hz_scaled = scale_vectors(Hx_top, Hz_top, arrow_length)
        
        ax.quiver(X, Z, Ex_scaled, Ez_scaled,
                color='blue', scale=1, scale_units='xy',
                width=0.006, headwidth=5)
        ax.quiver(X, Z, Hx_scaled, Hz_scaled,
                color='red', scale=1, scale_units='xy',
                width=0.006, headwidth=5)
    
        marker_size = 45
        ey_pos_mask = (Ey_top > threshold)
        ey_neg_mask = (Ey_top < -threshold)
        ax.scatter(X[ey_pos_mask], Z[ey_pos_mask], marker='o', s=marker_size, color='blue')
        ax.scatter(X[ey_neg_mask], Z[ey_neg_mask], marker='x', s=marker_size, color='blue')
        hy_pos_mask = (Hy_top > threshold)
        hy_neg_mask = (Hy_top < -threshold)
        ax.scatter(X[hy_pos_mask], Z[hy_pos_mask], marker='o', s=marker_size, color='red')
        ax.scatter(X[hy_neg_mask], Z[hy_neg_mask], marker='x', s=marker_size, color='red')

    def waveguide_side_section(self, mode, m, n, z, x_cut, ax):
        """Plot side (ZY-plane) field components at fixed X-cut position.
        
        Args:
            mode (int): 1=TE mode, 2=TM mode
            m (int): X-direction mode index
            n (int): Y-direction mode index
            z (float): Central Z-position in wavelengths (λ)
            x_cut (float): X-cut position (0-1 relative to width)
            ax (matplotlib.axes.Axes): Axes object for plotting
            
        Calculates:
            - Longitudinal field components (Ez, Hz)
            - Transverse components (Ey, Hy) using mode-specific relations
            - Includes both vector fields and scalar markers for:
              * Out-of-plane E-field components (blue markers)
              * Out-of-plane H-field components (red markers)
        """
        B = 2 * np.pi
        a, b = 1, 1
        z_vals = np.array([z - 0.5, z - 0.25, z, z + 0.25, z + 0.5])
        X_fixed = np.clip(x_cut * a, 0, a)
        regions_y = 2 * n + 1
        y_vals = np.linspace(0, b, regions_y)
        if n == 0:  # if n==0, center y at b/2
            y_vals = np.array([b/2])
        Y, Z = np.meshgrid(y_vals, z_vals, indexing='ij')

        siny = np.sin(n * np.pi * Y / b)
        cosy = np.cos(n * np.pi * Y / b)
        sinx_val = np.sin(m * np.pi * X_fixed / a)
        cosx_val = np.cos(m * np.pi * X_fixed / a)
        sinz = np.sin(B * Z)
        cosz = np.cos(B * Z)

        Ey_side, Ez_side = np.zeros_like(Y), np.zeros_like(Y)
        Hy_side, Hz_side = np.zeros_like(Y), np.zeros_like(Y)
        Ex_side, Hx_side = np.zeros_like(Y), np.zeros_like(Y)

        if mode == 1:  # TE Mode
            k_c = np.sqrt((m * np.pi / a)**2 + (n * np.pi / b)**2)
            scaling = B / k_c**2
            Hz_side = cosx_val * cosy * cosz
            Ex_side = -scaling * (n * np.pi / b) * cosx_val * siny * sinz
            Ey_side = -scaling * (m * np.pi / a) * sinx_val * cosy * sinz
            Hx_side = (m * np.pi / a) * sinx_val * cosy * sinz
            Hy_side = (n * np.pi / b) * cosx_val * siny * sinz

        elif mode == 2:  # TM Mode
            k_c = np.sqrt((m * np.pi / a)**2 + (n * np.pi / b)**2)
            scaling = B / k_c**2
            Ez_side = sinx_val * siny * cosz
            Ex_side = scaling * (m * np.pi / a) * cosx_val * siny * sinz
            Ey_side = -scaling * (n * np.pi / b) * sinx_val * cosy * sinz
            Hx_side = (n * np.pi / b) * sinx_val * cosy * sinz
            Hy_side = - (m * np.pi / a) * cosx_val * siny * sinz

        threshold = 1e-6
        Ey_side[np.abs(Ey_side) < threshold] = 0
        Hz_side[np.abs(Hz_side) < threshold] = 0
        Hy_side[np.abs(Hy_side) < threshold] = 0
        Hx_side[np.abs(Hx_side) < threshold] = 0

        ax.set_xticks(z_vals)
        ax.set_yticks(y_vals)
        ax.grid(True, linestyle='--', linewidth=0.5)
        ax.set_xlabel("Z-axis (propagation direction)")
        ax.set_ylabel("Y-axis")
        ax.set_title(f"Side Section (ZY) at X = {X_fixed:.2f}a")
        
        grid_spacing_y = b if regions_y == 1 else b / (regions_y - 1)  # avoid division by zero when n==0
        arrow_length = 0.4 * grid_spacing_y
        
        Ez_scaled, Ey_scaled = scale_vectors(Ez_side, Ey_side, arrow_length)
        Hz_scaled, Hy_scaled = scale_vectors(Hz_side, Hy_side, arrow_length)
        
        ax.quiver(Z, Y, Ez_scaled, Ey_scaled,
                color='blue', scale=1, scale_units='xy',
                width=0.006, headwidth=5)
        ax.quiver(Z, Y, Hz_scaled, Hy_scaled,
                color='red', scale=1, scale_units='xy',
                width=0.006, headwidth=5)
        
        marker_size = 60
        ex_pos_mask = (Ex_side > threshold)
        ex_neg_mask = (Ex_side < -threshold)
        ax.scatter(Z[ex_pos_mask], Y[ex_pos_mask], marker='o', s=marker_size, color='blue', edgecolor='black')
        ax.scatter(Z[ex_neg_mask], Y[ex_neg_mask], marker='x', s=marker_size, color='blue', linewidth=1.5)
        hx_pos_mask = (Hx_side < -threshold)
        hx_neg_mask = (Hx_side > threshold)
        ax.scatter(Z[hx_pos_mask], Y[hx_pos_mask], marker='o', s=marker_size, color='red', edgecolor='black')
        ax.scatter(Z[hx_neg_mask], Y[hx_neg_mask], marker='x', s=marker_size, color='red', linewidth=1.5)

if __name__ == "__main__":
    root = tk.Tk()
    root.state('zoomed')  # Open window fullscreen (maximized)
    app = WaveguideApp(root)
    root.protocol("WM_DELETE_WINDOW", app.on_close)
    root.mainloop()
