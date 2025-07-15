import matplotlib.pyplot as plt

def plot_multiplex_panels(panels, templates):
    fig, axes = plt.subplots(len(panels), 1, figsize=(8, 2*len(panels)))
    if len(panels)==1: axes=[axes]
    for ax, panel in zip(axes, panels):
        y = 10
        for p in panel:
            t = p['target']
            seq_len = len(templates[t])
            ax.set_xlim(0, seq_len)
            ax.set_ylim(0, 20)
            ax.set_yticks([10])
            ax.set_yticklabels([t])
            fs, fe = map(int, p['forward']['position'].split('-'))
            rs, re = map(int, p['reverse']['position'].split('-'))

            # plot forward as blue block
            ax.broken_barh([(fs, fe-fs)], (y, 5), facecolors='tab:blue')
            # plot reverse as green block
            ax.broken_barh([(rs, re-rs)], (y-6, 5), facecolors='tab:green')

        ax.grid(True)
    plt.tight_layout()
    plt.show()