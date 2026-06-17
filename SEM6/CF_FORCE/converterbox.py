import math

def convert_txt_to_xyz(input_file, output_file):

    with open(input_file) as f:
        lines = f.readlines()

    # ---- extract L ----
    L = None
    for line in lines:
        if line.startswith("# L="):
            L = float(line.split("=")[1])
            break

    if L is None:
        raise ValueError("Box size L not found in file")

    with open(output_file, "w") as f:

        i = 0
        while i < len(lines):

            if lines[i].startswith("FRAME"):
                i += 1
                particles = []

                while i < len(lines) and not lines[i].startswith("FRAME"):
                    if lines[i].strip() and not lines[i].startswith("#"):
                        particles.append(lines[i].strip())
                    i += 1

                N = len(particles)

                f.write(f"{N}\n")
                f.write(
                    f'Lattice="{L} 0 0  0 {L} 0  0 0 0.001" '
                    f'Origin="0 0 0" '
                    f'Properties=species:S:1:pos:R:3:theta:R:1:velo:R:3:Force:R:3\n'
                )

                for p in particles:
                    x, y, vx, vy, theta, fx, fy = map(float, p.split())

                    # write as 3D compatible XYZ (z = 0)
                    f.write(
                        f"A {x} {y} 0.0 "
                        f"{theta} "
                        f"{vx} {vy} 0.0 "
                        f"{fx} {fy} 0.0\n"
                    )

            else:
                i += 1


convert_txt_to_xyz("simulation_data.txt", "simulation.xyz")
print("Converted successfully")
