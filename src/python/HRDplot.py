import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt

#replace age_flame and mh_gspspec with your derived Bayesian ages and metallicities
#as well as the upper and lower bounds for uncertainties

def HRDage(csv_file, save_path, highlight_stars=None):
    df = pd.read_csv(csv_file, on_bad_lines='skip')
    df = df[df["lum_flame"] <= 5]

    parser = argparse.ArgumentParser(description="Plot HR Diagram with selectable y-axis.")
    parser.add_argument("--yaxis", choices=["lum_flame", "phot_g_mean_mag"], default="lum_flame",
                        help="Y-axis to plot: 'lum_flame' or 'phot_g_mean_mag'")
    args = parser.parse_known_args()[0]

    yaxis = args.yaxis
    ylabel = "Luminosity [L$_\odot$]" if yaxis == "lum_flame" else "G-band Magnitude"

    axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)[1]

    for ax in axes:
        ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in', top=True, right=True)
        ax.tick_params(axis='both', which='minor', length=4, width=1, direction='in', top=True, right=True)
        xticks_major = np.arange(3000, 7001, 500)
        xticks_major = xticks_major[1:-1]
        ax.set_xticks(xticks_major)
        ax.set_xticks(np.arange(3000, 7001, 100), minor=True)
        if yaxis == "lum_flame":
            yticks_major = np.arange(-0.5, 5.01, 0.5)
            yticks_major = yticks_major[1:-1]
            ax.set_yticks(yticks_major)
            ax.set_yticks(np.arange(-0.5, 5.01, 0.1), minor=True)
            # Change the font size of the actual numbers on the y axis here:
            ax.set_yticklabels([f"{tick:.1f}" for tick in yticks_major], minor=False, fontsize=10)
        else:
            yticks_major = np.arange(5, 16.1, 1)
            yticks_major = yticks_major[1:-1]
            ax.set_yticks(yticks_major)
            ax.set_yticks(np.arange(5, 16.1, 0.2), minor=True)

    sc1 = axes[0].scatter(
        df["teff_gspphot"], df[yaxis],
        c=df["age_flame"], cmap='viridis', s=20, edgecolor='k'
    )
    axes[0].set_xlabel("T$_{\\mathrm{eff}}$ [K]", fontsize=20)
    axes[0].set_xlim(3000, 7000)
    if yaxis == "lum_flame":
        axes[0].set_ylim(-0.5, 5)
    else:
        axes[0].set_ylim(16, 5)
    axes[0].set_ylabel(ylabel, fontsize=20)
    axes[0].invert_xaxis()
    cb1 = plt.colorbar(sc1, ax=axes[0], orientation='vertical', pad=0.02)
    cb1.ax.set_xlabel("Age Mean [Gyr]", labelpad=10, fontsize=14)
    cb1.ax.tick_params(labelsize=14)
    # Highlight specific stars with red boxes if provided
    if highlight_stars is not None and yaxis == "lum_flame":
        for star in highlight_stars:
            teff, lum = star
            box_width = 100   # K
            box_height = 0.1 # Lsun
            rect = plt.Rectangle(
                (teff - box_width/2, lum - box_height/2),
                box_width, box_height,
                linewidth=1.5, edgecolor='red', facecolor='none', zorder=10
            )
            axes[0].add_patch(rect)

    if yaxis == "lum_flame":
        box_x = [5500, 3000, 3000, 5500, 5500]
        box_y = [0, 0, 0.1, 0.1, 0]
        axes[1].plot(box_x, box_y, linestyle='--', color='black', alpha=0.5, linewidth=1)

    age_uncertainty = 0.5 * (df["age_flame_upper"] - df["age_flame_lower"])
    sc2 = axes[1].scatter(
        df["teff_gspphot"], df[yaxis],
        c=age_uncertainty, cmap='plasma', s=20, edgecolor='k'
    )
    axes[1].set_xlabel("T$_{\\mathrm{eff}}$ [K]", fontsize=20)
    axes[1].invert_xaxis()
    if yaxis == "phot_g_mean_mag":
        axes[1].set_ylim(16, 5)
    cb2 = plt.colorbar(sc2, ax=axes[1], orientation='vertical', pad=0.02)
    cb2.ax.set_xlabel("Age σ [Gyr]", labelpad=10, fontsize=14)
    cb2.ax.tick_params(labelsize=14)

    mask_cool_dim = (df["teff_gspphot"] < 5000) & (df["lum_flame"] < 0.1)
    avg_uncertainty_cool_dim = age_uncertainty[mask_cool_dim].mean()
    avg_uncertainty_rest = age_uncertainty[~mask_cool_dim].mean()

    #print(f"Average age uncertainty (teff < 5000K & L < 0.1): {avg_uncertainty_cool_dim:.3f} Gyr")
    #print(f"Average age uncertainty (rest): {avg_uncertainty_rest:.3f} Gyr")

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.show()

def HRDmet(csv_file, save_path):
    df = pd.read_csv(csv_file, on_bad_lines='skip')
    df = df[df["lum_flame"] <= 5]
    df = df[df["age_flame"].notna()]

    parser = argparse.ArgumentParser(description="Plot HR Diagram with selectable y-axis.")
    parser.add_argument("--yaxis", choices=["lum_flame", "phot_g_mean_mag"], default="lum_flame",
                        help="Y-axis to plot: 'lum_flame' or 'phot_g_mean_mag'")
    args, unknown = parser.parse_known_args()

    yaxis = args.yaxis
    ylabel = "Luminosity [L$_\odot$]" if yaxis == "lum_flame" else "G-band Magnitude"

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for ax in axes:
        ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in', top=True, right=True)
        ax.tick_params(axis='both', which='minor', length=4, width=1, direction='in', top=True, right=True)
        xticks_major = np.arange(3000, 7001, 500)
        xticks_major = xticks_major[1:-1]
        ax.set_xticks(xticks_major)
        ax.set_xticks(np.arange(3000, 7001, 100), minor=True)
        if yaxis == "lum_flame":
            yticks_major = np.arange(-0.5, 5.01, 0.5)
            yticks_major = yticks_major[1:-1]
            ax.set_yticks(yticks_major)
            ax.set_yticks(np.arange(-0.5, 5.01, 0.1), minor=True)
            ax.set_yticklabels([f"{tick:.1f}" for tick in yticks_major], minor=False, fontsize=10)
        else:
            yticks_major = np.arange(5, 16.1, 1)
            yticks_major = yticks_major[1:-1]
            ax.set_yticks(yticks_major)
            ax.set_yticks(np.arange(5, 16.1, 0.2), minor=True)

    sc1 = axes[0].scatter(
        df["teff_gspphot"], df[yaxis],
        c=df["mh_gspspec"], cmap='viridis', s=20, edgecolor='k', vmin=-2
    )
    axes[0].set_xlabel("T$_{\\mathrm{eff}}$ [K]", fontsize=20)
    axes[0].set_xlim(3000, 7000)
    if yaxis == "lum_flame":
        axes[0].set_ylim(-0.5, 5)
    else:
        axes[0].set_ylim(16, 5)
    axes[0].set_ylabel(ylabel, fontsize=20)
    axes[0].invert_xaxis()
    cb1 = plt.colorbar(sc1, ax=axes[0], orientation='vertical', pad=0.02)
    cb1.ax.set_xlabel("[M/H] Mean", labelpad=10, fontsize=14)
    cb1.ax.tick_params(labelsize=14)

    met_uncertainty = 0.5 * (df["mh_gspspec_upper"] - df["mh_gspspec_lower"])
    sc2 = axes[1].scatter(
        df["teff_gspphot"], df[yaxis],
        c=met_uncertainty, cmap='plasma', s=20, edgecolor='k'
    )
    axes[1].set_xlabel("T$_{\\mathrm{eff}}$ [K]", fontsize=20)
    axes[1].invert_xaxis()
    if yaxis == "phot_g_mean_mag":
        axes[1].set_ylim(16, 5)
    cb2 = plt.colorbar(sc2, ax=axes[1], orientation='vertical', pad=0.02)
    cb2.ax.set_xlabel("[M/H] σ", labelpad=10, fontsize=14)
    cb2.ax.tick_params(labelsize=14)

    mask_cool_dim = (df["teff_gspphot"] < 5000) & (df["lum_flame"] < 0.1)
    avg_uncertainty_cool_dim = met_uncertainty[mask_cool_dim].mean()
    avg_uncertainty_rest = met_uncertainty[~mask_cool_dim].mean()

    print(f"Average metallicity uncertainty (teff < 5000K & L < 0.1): {avg_uncertainty_cool_dim:.3f}")
    print(f"Average metallicity uncertainty (rest): {avg_uncertainty_rest:.3f}")

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.show()


if __name__ == "__main__":
    highlight_stars = [(4607, 0.202), (6465, 4.38), (5896, 1.356), (5536,0.745) ] #(6000, 2.5), (3500, 0.05)]
    csv_file = "/path/to/your/csv/file"
    #HRDage(csv_file, "/path/to/your/savefigfile/HRD_age_colored.png", highlight_stars=highlight_stars)
    HRDmet(csv_file, "/path/to/your/savefigfile//HRD_met_colored.png")
    
