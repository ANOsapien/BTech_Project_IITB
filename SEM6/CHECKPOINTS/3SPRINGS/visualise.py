#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import hsv_to_rgb
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.animation import FFMpegWriter


def read_simulation_data(filename):
    """Read simulation data from file"""
    frames = []
    current_frame = None
    L = None
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                if line.startswith('# N='):
                    parts = line.split()
                    for p in parts:
                        if p.startswith('L='):
                            L = float(p.split('=')[1])
                continue
            
            if line.startswith('FRAME'):
                if current_frame is not None:
                    frames.append(current_frame)
                current_frame = {'particles': []}
            else:
                data = list(map(float, line.split()))
                current_frame['particles'].append({
                    'x': data[0],
                    'y': data[1],
                    'theta': data[2],
                    'vx': data[3],
                    'vy': data[4]
                })
        
        if current_frame is not None:
            frames.append(current_frame)
    
    return frames, L

def theta_to_color(theta):
    """Convert angle theta to RGB color using HSV"""
    h = (theta % (2 * np.pi)) / (2 * np.pi)
    s = 0.9
    v = 1.0
    return hsv_to_rgb([h, s, v])

def setup_figure(L):
    """Create figure and axis with colorbar"""
    fig, ax = plt.subplots(figsize=(11, 10))
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=11)
    ax.set_ylabel('y', fontsize=11)
    
    # Create colorbar
    norm = Normalize(vmin=0, vmax=2*np.pi)
    sm = ScalarMappable(cmap='hsv', norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label='Direction (radians)', pad=0.02)
    cbar.set_ticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    cbar.set_ticklabels(['0', 'π/2', 'π', '3π/2', '2π'])
    
    return fig, ax

def update_frame(frame_data, ax, L, arrow_scale=0.25):
    """Update plot for a frame"""
    ax.clear()
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=11)
    ax.set_ylabel('y', fontsize=11)
    ax.set_title(f'Frame {frame_data["frame_num"]}', fontsize=12, fontweight='bold')
    
    particles = frame_data['particles']
    
    # Plot all particles and arrows
    for p in particles:
        color = theta_to_color(p['theta'])
        
        # Plot particle as filled circle (no edge for smoothness)
        ax.plot(p['x'], p['y'], 'o', color=color, markersize=6.5, markeredgewidth=0)
        
        # Draw direction arrow with matching color
        ax.arrow(p['x'], p['y'], 
                arrow_len := arrow_scale * np.cos(p['theta']), 
                arrow_scale * np.sin(p['theta']),
                head_width=0.12, head_length=0.1, 
                fc=color, ec=color, alpha=0.8, linewidth=1.2)

def animate_simulation(filename, L):
    """Create animation from simulation data"""
    frames, L_from_file = read_simulation_data(filename)
    
    if L is None:
        L = L_from_file
    
    if not frames:
        print("No frames found in data file!")
        return
    
    print(f"Loaded {len(frames)} frames, L={L}")
    
    fig, ax = setup_figure(L)
    
    # Add frame information to each frame
    for i, frame in enumerate(frames):
        frame['frame_num'] = i
    
    def update_anim(frame_idx):
        update_frame(frames[frame_idx], ax, L, arrow_scale=0.22)
        return ax,
    
    # Smoother animation: faster frame rate, shorter interval
    anim = animation.FuncAnimation(fig, update_anim, frames=len(frames),
                                  interval=10, repeat=True, blit=False)
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    filename = 'simulation_data.txt'
    
    # ANIMATION SPEED CONTROL
    # Decrease for faster animation, increase for slower
    ANIMATION_INTERVAL_MS = 0.5  # milliseconds between frames
    
    frames, L = read_simulation_data(filename)
    
    if not frames:
        print("No frames found in data file!")
        exit(1)
    
    print(f"Loaded {len(frames)} frames, L={L}")
    print(f"Total duration: {len(frames) * ANIMATION_INTERVAL_MS / 1000:.1f} seconds")
    
    fig, ax = setup_figure(L)
    
    # Add frame information
    for i, frame in enumerate(frames):
        frame['frame_num'] = i
    
    def update_anim(frame_idx):
        update_frame(frames[frame_idx], ax, L, arrow_scale=0.22)
        return ax,
    
    anim = animation.FuncAnimation(fig, update_anim, frames=len(frames),
                                  interval=ANIMATION_INTERVAL_MS, repeat=True, blit=False)


# Save animation as MP4 with faster playback speed
    writer = FFMpegWriter(fps=30)  # try fps=60 or even 120
    anim.save("simulation_fast.mp4", writer=writer)
    print("Saved high-speed video as simulation_fast.mp4")

    
    plt.tight_layout()
    plt.show()