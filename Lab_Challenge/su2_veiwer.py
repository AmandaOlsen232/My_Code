def read_su2(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    nodes = []
    elements = []
    boundaries = {}

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith('%'):
            i += 1
            continue

        if line.startswith('NPOIN='):
            n_nodes = int(line.split('=')[1])
            i += 1
            while i < len(lines):
                line = lines[i].strip()
                if not line or line.startswith('%'):
                    i += 1
                    continue
                if line.startswith('NELEM=') or line.startswith('NMARK='): # If we reach the next section, break
                    break
                parts = line.split()
                try:
                    x, y = float(parts[0]), float(parts[1])
                    nodes.append((x, y))
                except ValueError:
                    raise ValueError(f"Invalid float conversion in node line: '{line}'")
                i += 1
            continue

        elif line.startswith('NELEM='):
            n_elements = int(line.split('=')[1])
            i += 1
            count = 0
            while count < n_elements and i < len(lines):
                line = lines[i].strip()
                if not line or line.startswith('%'):
                    i += 1
                    continue
                parts = list(map(int, line.strip().split()))
                etype = parts[0]
                if etype == 5:
                    elements.append(parts[1:4])
                elif etype == 9:
                    elements.append(parts[1:5])
                else:
                    print(f"Unsupported element type {etype} at line {i}")
                count += 1
                i += 1
            continue

        elif line.startswith('NMARK='):
            n_mark = int(line.split('=')[1])
            i += 1
            for _ in range(n_mark):
                while lines[i].strip().startswith('%') or not lines[i].strip():
                    i += 1
                tag = lines[i].strip().split('=')[1].strip()
                i += 1
                while lines[i].strip().startswith('%') or not lines[i].strip():
                    i += 1
                n_edges = int(lines[i].strip().split('=')[1])
                i += 1
                edges = []
                for _ in range(n_edges):
                    while lines[i].strip().startswith('%') or not lines[i].strip():
                        i += 1
                    parts = list(map(int, lines[i].strip().split()))
                    if parts[0] == 3:
                        edges.append((parts[1], parts[2]))
                    i += 1
                boundaries[tag] = edges
            continue

        i += 1

    return nodes, elements, boundaries

def plot_su2_mesh(nodes, elements, boundaries):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_facecolor('#F5F5DC')
    for elem in elements:
        poly = Polygon([nodes[i] for i in elem], closed=True, edgecolor='#8B5C2D', facecolor='#A3C585', alpha=0.7)
        ax.add_patch(poly)

    xs, ys = zip(*nodes)
    ax.plot(xs, ys, 'o', markersize=1, color='#4F97A3', alpha=0.5)

    # Plot boundaries in different colors
    boundary_colors = ['#D7263D', '#1CA9C9', '#FFD700', '#6B8E23', '#FF8C00']
    for idx, (tag, edges) in enumerate(boundaries.items()):
        color = boundary_colors[idx % len(boundary_colors)]
        for n0, n1 in edges:
            x_vals = [nodes[n0][0], nodes[n1][0]]
            y_vals = [nodes[n0][1], nodes[n1][1]]
            ax.plot(x_vals, y_vals, color=color, linewidth=2, label=tag if n0 == edges[0][0] and n1 == edges[0][1] else "")

    ax.set_aspect('equal')
    ax.set_title('SU2 Mesh Visualization')
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())
    plt.grid(True, color='#FFD700', alpha=0.2)
    plt.show()

# Example usage
nodes, elements, boundaries = read_su2("mesh_NACA2412_inv.su2")
plot_su2_mesh(nodes, elements, boundaries)
