import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from tqdm import tqdm
import re

def natural_sort_key(file_name):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', file_name)]

def animate_structure(nodesCoordFile, edgesFile, dataFolder, outputFile):
    # Read nodes
    nodes = []
    with open(nodesCoordFile, 'r') as f:
        for line in f:
            x, y = map(float, line.strip().split(','))
            nodes.append((x, y))

    # Read edges
    edges = []
    with open(edgesFile, 'r') as f:
        for line in f:
            idx1, idx2 = map(int, line.strip().split(','))
            edges.append((idx1 - 1, idx2 - 1))  # Convert to 0-based indexing

    # Get all deformation files
    deformation_files = sorted([os.path.join(dataFolder, file) for file in os.listdir(dataFolder) if file.startswith("time") and file.endswith(".txt") and int(file[4:-4]) % 100 == 0], key=natural_sort_key)
    
    print(deformation_files)
    
    print(f"Found {len(deformation_files)} deformation files.")
    # Initialize the plot
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim(min(x for x, y in nodes) - 0.1, max(x for x, y in nodes) + 0.1)
    ax.set_ylim(min(y for x, y in nodes) - 0.1, max(y for x, y in nodes) + 0.1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Animation de la structure (Facteur d\'amplification = 300)')
    ax.grid(True)
    ax.axis('equal')

    lines = []
    for _ in edges:
        line, = ax.plot([], [], color='blue', linewidth=0.8)
        lines.append(line)

    def update(frame):
        # Read displacements for the current frame
        displacements = []
        amplification_factor = 300
        with open(deformation_files[frame], 'r') as f:
            for line in f:
                try :
                    ux, uy, _, _ = map(float, line.strip().split())
                    displacements.append((ux * amplification_factor, uy * amplification_factor))
                except :
                    print(f"Error reading line in {deformation_files[frame]}: {line.strip()}")
                    exit(1)
        # Update edges with displacements
        for i, (idx1, idx2) in enumerate(edges):
            x_values = [nodes[idx1][0] + displacements[idx1][0], nodes[idx2][0] + displacements[idx2][0]]
            y_values = [nodes[idx1][1] + displacements[idx1][1], nodes[idx2][1] + displacements[idx2][1]]
            lines[i].set_data(x_values, y_values)

        return lines

    ani = animation.FuncAnimation(fig, update, frames=len(deformation_files), blit=True, repeat=False)

    # Save the animation as a GIF instead of MP4
    ani.save(outputFile.replace('.mp4', '.gif'), writer='pillow', fps=30)
    plt.close()

# Example usage
animate_structure('./modelCoordinates.txt', './edges.txt', './plots/data', 'structure_animation.mp4')
