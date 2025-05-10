import matplotlib.pyplot as plt

def plot_structure(nodesCoordFile, edgesFile, displacement, amplificationFactor):
    # On recup les noeuds
    nodes = []
    with open(nodesCoordFile, 'r') as f:
        for line in f:
            x, y = map(float, line.strip().split(','))
            nodes.append((x, y))

    # Read displacements from final.txt
    displacements = []
    with open(displacement, 'r') as f:
        for line in f:
            ux, uy, _, _ = map(float, line.strip().split())
            displacements.append((ux * amplificationFactor, uy * amplificationFactor))
    # On recup le contour
    edges = []
    with open(edgesFile, 'r') as f:
        for line in f:
            idx1, idx2 = map(int, line.strip().split(','))
            edges.append((idx1 - 1, idx2 - 1))

    plt.figure(figsize=(10, 10))

    # Plot edges with displacements
    for idx1, idx2 in edges:
        x_values = [nodes[idx1][0] + displacements[idx1][0], nodes[idx2][0] + displacements[idx2][0]]
        y_values = [nodes[idx1][1] + displacements[idx1][1], nodes[idx2][1] + displacements[idx2][1]]
        plt.plot(x_values, y_values, color='green', linewidth=0.8, linestyle='--', label='Displaced Structure' if idx1 == 0 else "")

    # Plot original structure edges
    for idx1, idx2 in edges:
        x_values = [nodes[idx1][0], nodes[idx2][0]]
        y_values = [nodes[idx1][1], nodes[idx2][1]]
        plt.plot(x_values, y_values, color='blue', linewidth=0.8, label='Original Structure' if idx1 == 0 else "")

    # Plot nodes
    for x, y in nodes:
        plt.scatter(x, y, color='red', s=1, label='Original Nodes' if nodes.index((x, y)) == 0 else "")

    # Plot displaced nodes
    for i, (x, y) in enumerate(nodes):
        plt.scatter(x + displacements[i][0], y + displacements[i][1], color='orange', s=10, label='Displaced Nodes' if i == 0 else "")

    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title('Structure Plot with Displacement')
    plt.grid(True)
    plt.axis('equal')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    nodesCoordFile = './modelCoordinates.txt'
    edgesFile = './edges.txt'
    amplificationFactor = 300
    plot_structure(nodesCoordFile, edgesFile, './final.txt', amplificationFactor)
    