import numpy as np

import matplotlib.pyplot as plt

def load_posterior_from_text(filepath):
    # Assumes the file contains a flat list of probability values
    data = np.loadtxt(filepath)
    # Reshape to 50 (age) x 48 (metallicity) grid
    return data.reshape((50, 48))

def plot_posterior_age_metallicity(posterior, min_age, max_age, min_mh, max_mh, savelocation, logTe, logL, BV, Mv, logscale=True):
    plt.figure(figsize=(8, 6))
    extent = [min_mh, max_mh, min_age, max_age]
    plot_data = posterior.copy()
    if logscale:
        min_nonzero = plot_data[plot_data > 0].min() if np.any(plot_data > 0) else 1e-10
        fill_value = min_nonzero / 10.0
        plot_data[plot_data == 0] = fill_value
        norm = plt.matplotlib.colors.LogNorm()
    else:
        norm = None
    plt.imshow(plot_data, extent=extent, origin='lower', aspect='auto', cmap='viridis', norm=norm)
    plt.colorbar(label='Posterior Probability')
    plt.xlabel('Metallicity [M/H]')
    plt.ylabel('Age [Gyr]')
    legend_text = f"logTe = {logTe:.2f}\nlogL = {logL:.2f}\nGbp-Grp = {BV:.2f}\nGmag = {Mv:.2f}"
    x_text = min_mh + 0.05 * (max_mh - min_mh)
    y_text = max_age - 0.05 * (max_age - min_age)
    plt.text(
        x=x_text,
        y=y_text,
        s=legend_text,
        fontsize=10,
        verticalalignment='top',
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='black')
    )
    plt.savefig(savelocation)
    plt.show()
    
def plot_4_posteriors_grid(param_list, basepath, min_age, max_age, min_mh, max_mh, logscale=True, savelocation=None):
    """
    param_list: list of dicts, each with keys 'logTe', 'logL', 'BV', 'Mv'
    basepath: directory where posterior files are stored
    min_age, max_age, min_mh, max_mh: axis limits
    logscale: use log scale for color
    savelocation: path to save the figure (optional)
    """
    # Article/report scale: slightly smaller, more compact
    fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True, sharey=True)
    extent = [min_mh, max_mh, min_age, max_age]
    ims = []
    for idx, params in enumerate(param_list):
        row, col = divmod(idx, 2)
        ax = axes[row, col]
        filename = (
            f"posterior_manual_logTe{params['logTe']:.2f}_logLLo{params['logL']:.2f}_GbpGrp{params['BV']:.2f}_G{params['Mv']:.2f}.txt"
        )
        filepath = basepath + filename
        posterior = load_posterior_from_text(filepath)
        plot_data = posterior.copy()
        if logscale:
            min_nonzero = plot_data[plot_data > 0].min() if np.any(plot_data > 0) else 1e-10
            fill_value = min_nonzero / 10.0
            plot_data[plot_data == 0] = fill_value
            norm = plt.matplotlib.colors.LogNorm()
        else:
            norm = None
        im = ax.imshow(plot_data, extent=extent, origin='lower', aspect='auto', cmap='viridis', norm=norm)
        ims.append(im)
        legend_text = (
            f"$\\log T_{{\\mathrm{{eff}}}}$ = {params['logTe']:.2f}\n"
            f"$\\log (L/L_\\odot)$ = {params['logL']:.2f}\n"
            f"$G_{{\\mathrm{{BP}}}} - G_{{\\mathrm{{RP}}}}$ = {params['BV']:.2f}\n"
            f"$G$ = {params['Mv']:.2f}"
        )
        x_text = min_mh + 0.05 * (max_mh - min_mh)
        y_text = max_age - 0.05 * (max_age - min_age)
        ax.text(
            x_text, y_text, legend_text, fontsize=11, verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='black')
        )
        if row == 1:
            ax.set_xlabel('Metallicity [M/H]', fontsize=18)
        if col == 0:
            ax.set_ylabel('Age [Gyr]', fontsize=18)
    # Place colorbar on the right outside the grid, label below
    fig.tight_layout(rect=[0, 0, 0.92, 1])
    cbar = fig.colorbar(ims[0], ax=axes, orientation='vertical', fraction=0.04, pad=0.02)
    cbar.ax.set_xlabel('P(Age, [M/H])', labelpad=10, fontsize=11)
    if savelocation:
        plt.savefig(savelocation, bbox_inches='tight', dpi=300)
    plt.show()

# Example: define the 4 sets of parameters here
param_list = [
    {'logTe': 3.66, 'logL': -0.69, 'BV': 1.31, 'Mv': 10.48},
    {'logTe': 3.74, 'logL': -0.13, 'BV': 0.93, 'Mv': 9.55},
    {'logTe': 3.77, 'logL': 0.13, 'BV': 0.74, 'Mv': 8.98},
    {'logTe': 3.81, 'logL': 0.64, 'BV': 0.58, 'Mv': 8.06},
]
if __name__ == "__main__":
    # Define parameters
    logTe = 3.74
    logL = -0.13
    BV = 0.93
    Mv = 9.55

    # File path construction
    basepath = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/outputdata/posteriors/"
    filename = f"posterior_manual_logTe{logTe:.2f}_logLLo{logL:.2f}_GbpGrp{BV:.2f}_G{Mv:.2f}.txt"
    filepath = basepath + filename

    # Load posterior
    posterior = load_posterior_from_text(filepath)

    # Plot settings
    min_age = 0.1
    max_age = 14.0
    min_mh = -1.5
    max_mh = 0.42
    logscale = True
    logscale_str = "_logscale" if logscale else ""
    savelocation = (
        "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/plots/posteriors/"
        f"posterior{logTe:.2f}_{logL:.2f}_{BV:.2f}_{Mv:.2f}{logscale_str}.png"
    )

    #plot_posterior_age_metallicity(
    #    posterior, min_age, max_age, min_mh, max_mh,
    #    savelocation, logTe, logL, BV, Mv, logscale=logscale
    #)
    
    plot_4_posteriors_grid(
        param_list, basepath, min_age, max_age, min_mh, max_mh,
        logscale=logscale,
        savelocation="/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/plots/posteriors/posterior_grid.png"
    )
    