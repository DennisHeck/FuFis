import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import networkx as nx
from matplotlib import cm, colors, colorbar


def directed_network(edge_df, edge_cols, edge_colour_col=None, node_values=None, node_cmap='bwr',
                     node_cmap_centre=False, node_cbar_label='', xsize=12, ysize=10, title='', plot_out=''):

    all_edges = edge_df[edge_cols].values
    G = nx.DiGraph()  # create empty graph
    G.add_edges_from(all_edges)

    if edge_colour_col:
        edge_colour_dict = {e[0]+"#"+e[1]: e[2] for e in edge_df[edge_cols+[edge_colour_col]].values}
        edge_colours = [edge_colour_dict['#'.join(e)] for e in G.edges]

    if node_values:
        cmap = cm.get_cmap(node_cmap)
        if node_cmap_centre or node_cmap_centre == 0:
            boundup = max([abs(min(node_values.values())), max(node_values.values())])
            boundlow = -boundup
        else:
            boundup = max(node_values.values())
            boundlow = min(node_values.values())
        norm = plt.Normalize(boundlow, boundup)
        node_colours = [to_hex(cmap(norm(node_values[node]))) for node in G.nodes]

    f, ax = plt.subplots(figsize=(xsize, ysize))
    pos = nx.circular_layout(G)
    nx.draw_networkx(G, connectionstyle='arc3, rad = 0.03', with_labels=True, arrows=True,
                     edge_color=None if not edge_colour_col else edge_colours,
                     node_color=None if not node_values else node_colours,
                     pos=pos, width=0.5, node_size=600, arrowsize=9, font_size=6)
    for e, edge in enumerate(G.edges):  # Manually add self-loops.
        if edge[0] == edge[1]:
            ax.annotate("",
                        xy=[pos[edge[0]][0], pos[edge[0]][1 ] -0.055], xycoords='data',
                        xytext=[pos[edge[0]][0], pos[edge[0]][1 ] +0.055], textcoords='data',
                        arrowprops=dict(arrowstyle="-|>", color=edge_colours[e], linewidth=0.5,
                                        shrinkA=0, shrinkB=0,
                                        patchA=None, patchB=None,
                                        connectionstyle="arc3, rad=1.5",
                                        ),
                        )
    plt.title(title)
    if node_values:
        # Add a manual colormap to be able to edit it.
        sm = plt.cm.ScalarMappable(cmap=cm.get_cmap(node_cmap), norm=colors.Normalize(vmin=boundlow, vmax=boundup))
        plt.colorbar(sm)
        sm.colorbar.set_label(node_cbar_label)
    plt.savefig(plot_out + "_DirNetwork.pdf")
    plt.close()

