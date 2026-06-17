import os

# Root folder containing all simulation folders
root_dir = "results"

import re

import re

def convert_txt_to_xyz(input_file, output_file):
    """
    Convert simulation data from custom TXT format to XYZ format for OVITO.
    
    XYZ format:
    Line 1: Number of atoms
    Line 2: Comment line (can contain parameters)
    Lines 3+: atom_type x y z vx vy vz (or just x y z)
    """
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Parse parameters from first line
    params_line = lines[0].strip()
    params = {}
    for param in params_line.replace('# Simulation Parameters', '').split():
        if '=' in param:
            key, val = param.split('=')
            params[key] = val
    
    N = int(params.get('N', 0))
    
    with open(output_file, 'w') as f:
        frame_num = 0
        i = 1
        
        while i < len(lines):
            line = lines[i].strip()
            
            # Check if it's a FRAME line
            if line.startswith('FRAME'):
                frame_num = int(line.split()[1])
                
                # Count particles in this frame
                particle_lines = []
                i += 1
                while i < len(lines) and not lines[i].strip().startswith('FRAME'):
                    if lines[i].strip():  # Skip empty lines
                        particle_lines.append(lines[i].strip())
                    i += 1
                
                # Write XYZ format
                f.write(f"{len(particle_lines)}\n")
                f.write(f"Frame {frame_num}, {params_line}\n")
                
                # Write particle data (format: x y theta vx vy)
                # Convert to 3D XYZ: x y z=0 theta vx vy
                for idx, p_line in enumerate(particle_lines):
                    values = p_line.split()
                    if len(values) != 5:
                        print(f"Warning: Frame {frame_num}, particle {idx+1} has {len(values)} values: {p_line}")
                        continue  # Skip malformed lines
                    x, y, theta, vx, vy = values
                    # Format: type x y z theta vx vy
                    # Set z=0 for 2D visualization, keep theta as extra column
                    f.write(f"1 {x} {y} 0.0 {theta} {vx} {vy}\n")
            else:
                i += 1




input_file =  "simulation_data.txt"
output_file =  "simulation.xyz"

convert_txt_to_xyz(input_file, output_file)
