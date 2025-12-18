import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import medfilt
import matplotlib.pyplot as plt



# === Configuration ===
csv_file = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/outputdata/velocitydispersion/sigma_agemean_metmean_binned.csv"
csv_file2 = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/outputdata/velocitydispersion/combined_sigma_age_binned_age_verror_parallax_MScuts100.csv"
csv_file3 = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/outputdata/velocitydispersion/combined_sigma_age_binned_agecut_verrorcut100.csv"
mhcsv_file = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/outputdata/velocitydispersion/sigma_mh_binned_agecut_maxage100.csv"
filter_zeros = True     # ignore bins where all dispersions are 0
plot_all = True         # if False, only plot one component below
single_component = "sigma_total" # other options: sigma_U, sigma_V, sigma_W
cmap = "viridis"
parameter_names = "teff_luminosity"  # used for saving the plot
filepath = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/"
components = ["sigma_U", "sigma_V", "sigma_W", "sigma_total"]
labels = {
    "sigma_U": r"$\sigma_\mathrm{U}$ [km/s]",
    "sigma_V": r"$\sigma_\mathrm{V}$ [km/s]",
    "sigma_W": r"$\sigma_\mathrm{W}$ [km/s]",
    "sigma_total": r"$\sigma_\mathrm{total}$ [km/s]",
    "title": r"$\sigma$ Components for binned $T_\mathrm{eff}$ and Luminosity"
}

def load_data(csv_file_path, filter_zeros=True, components=None):
    df = pd.read_csv(csv_file_path)
    if filter_zeros and components is not None:
        df = df[(df[components] != 0).any(axis=1)]
    return df

def plot_all_components(df, components, labels, cmap, filepath, parameter_names):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    vmin = min(df[comp].min() for comp in components)
    vmax = max(df[comp].max() for comp in components)
    im = None
    for ax, comp in zip(axes, components):
        pivot = df.pivot_table(
            index="param2_center", columns="param1_center", values=comp, aggfunc="first"
        )
        x = pivot.columns.values
        y = pivot.index.values
        X, Y = np.meshgrid(x, y)
        Z = pivot.values
        im = ax.pcolormesh(X, Y, Z, shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(labels[comp])
        # Only set x-axis label for bottom row plots (indices 2 and 3)
        if ax in [axes[2], axes[3]]:
            ax.set_xlabel("Effective Temperature [K]")
        else:
            ax.set_xlabel("")
        # Only set y-axis label for left column plots (indices 0 and 2)
        if ax in [axes[0], axes[2]]:
            ax.set_ylabel("Luminosity [L☉]")
        else:
            ax.set_ylabel("")
    fig.colorbar(im, ax=axes, orientation='vertical', fraction=0.025, pad=0.04)
    plt.suptitle(labels["title"], fontsize=16)
    plt.tight_layout()
    plt.savefig(filepath + f"/plots/dispersion/velocity_dispersion_{parameter_names}.png", dpi=300)
    plt.show()

def plot_single_component(df, component, labels, cmap, filepath):
    if component not in labels:
        raise ValueError(f"Invalid component name: {component}")
    pivot = df.pivot_table(
        index="param2_center", columns="param1_center", values=component, aggfunc="first"
    )
    x = pivot.columns.values
    y = pivot.index.values
    Z = pivot.values
    X, Y = np.meshgrid(x, y)
    plt.figure(figsize=(8, 6))
    im = plt.pcolormesh(X, Y, Z, shading='auto', cmap=cmap)
    plt.title(labels[component])
    plt.xlabel("Age [Gyr]")
    plt.ylabel("[M/H]")
    plt.colorbar(im)
    plt.tight_layout()
    plt.savefig(filepath + f"/plots/dispersion/velocity_dispersion_{component}.png", dpi=300)
    plt.show()

def plot_dispersion_vs_age(df, filepath, plot_loglog=False):
    age_col = "param_center"
    disp_col = "sigma_total"
    err_col = "error_sigma_total"
    # Group by age and compute mean velocity dispersion and error for each age bin
    age_disp = df.groupby(age_col)[[disp_col, err_col]].mean().reset_index()
    # Manually specify the bin indices (row numbers) to exclude
    bins_to_exclude = list(range(0, 8)) + list(range(72, 74)) + [60, 81, 88, 96]
    mask = ~age_disp.index.isin(bins_to_exclude)
    age_disp_filtered = age_disp[mask].reset_index(drop=True)
    # Define a power law function: sigma = sigma_0 * age^beta
    def power_law(age, sigma_0, beta):
        return sigma_0 * np.power(age, beta)
    # Use the errors as sigma in curve_fit to weight the fit and get uncertainties
    popt, pcov = curve_fit(
        power_law,
        age_disp_filtered[age_col],
        age_disp_filtered[disp_col],
        sigma=age_disp_filtered[err_col],
        absolute_sigma=True,
        p0=[1, 0.5]
    )
    sigma_0_fit, beta_fit = popt
    perr = np.sqrt(np.diag(pcov))
    sigma_0_err, beta_err = perr
    plt.figure(figsize=(7, 5))
    scatter_color = "#1f77b4"
    fit_color = "#d62728"
    age_fit = np.linspace(age_disp_filtered[age_col].min(), age_disp_filtered[age_col].max(), 200)
    if plot_loglog:
        plt.errorbar(
            age_disp_filtered[age_col], age_disp_filtered[disp_col],
            yerr=age_disp_filtered[err_col],
            fmt='o', color=scatter_color, ecolor='gray', elinewidth=1.5, capsize=3,
            label=r"Mean $\sigma_\mathrm{total}$ per Age bin"
        )
        plt.plot(
            age_fit, power_law(age_fit, *popt),
            color=fit_color, linestyle=':', linewidth=2,
            label=rf"Fit: $\sigma = {sigma_0_fit:.2f} \pm {sigma_0_err:.2f} \times \tau^{{{beta_fit:.2f} \pm {beta_err:.2f}}}$"
        )
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("Age [Gyr] (log)", fontsize=13)
        plt.ylabel(r"$\sigma_\mathrm{total}$ [km/s] (log)", fontsize=13)
        plt.title(r"Velocity Dispersion vs Age (log-log)", fontsize=15)
    else:
        plt.errorbar(
            age_disp_filtered[age_col], age_disp_filtered[disp_col],
            yerr=age_disp_filtered[err_col],
            fmt='o', color=scatter_color, ecolor='gray', elinewidth=1.5, capsize=3,
            label=r"Mean $\sigma_\mathrm{total}$ per Age bin"
        )
        plt.plot(
            age_fit, power_law(age_fit, *popt),
            color=fit_color, linestyle=':', linewidth=2,
            label=rf"Fit: $\sigma = {sigma_0_fit:.2f} \pm {sigma_0_err:.2f} \times \tau^{{{beta_fit:.2f} \pm {beta_err:.2f}}}$"
        )
        plt.xlabel("Age [Gyr]", fontsize=13)
        plt.ylabel(r"$\sigma_\mathrm{total}$ [km/s]", fontsize=13)
        plt.title(r"Velocity Dispersion vs Age", fontsize=15)
    plt.grid(True, which='both', linestyle='--', alpha=0.3)
    plt.legend(fontsize=12, frameon=True)
    plt.tight_layout()
    plt.savefig(filepath + f"/plots/dispersion/velocity_dispersion_vs_age_powerlaw100bins{'_loglog' if plot_loglog else ''}.png", dpi=300)
    plt.show()

def plot_dispersion_vs_age_all_components(df, filepath, plot_loglog=False):

    age_col = "param_center"
    err_cols = {
        "sigma_U": "error_sigma_U",
        "sigma_V": "error_sigma_V",
        "sigma_W": "error_sigma_W",
        "sigma_total": "error_sigma_total"
    }
    components = ["sigma_U", "sigma_V", "sigma_W", "sigma_total"]
    labels = {
        "sigma_U": r"$\sigma_\mathrm{U}$ [km/s]",
        "sigma_V": r"$\sigma_\mathrm{V}$ [km/s]",
        "sigma_W": r"$\sigma_\mathrm{W}$ [km/s]",
        "sigma_total": r"$\sigma_\mathrm{total}$ [km/s]"
    }
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    bins_to_exclude = list(range(0, 8)) + list(range(72, 74)) + [60, 81, 88, 96]
    scatter_color = "#1f77b4"
    fit_color = "#d62728"

    def power_law(age, sigma_0, beta):
        return sigma_0 * np.power(age, beta)

    for i, comp in enumerate(components):
        disp_col = comp
        err_col = err_cols[comp]
        age_disp = df.groupby(age_col)[[disp_col, err_col]].mean().reset_index()
        mask = ~age_disp.index.isin(bins_to_exclude)
        age_disp_filtered = age_disp[mask].reset_index(drop=True)
        # Use the errors as sigma in curve_fit to weight the fit and get uncertainties
        popt, pcov = curve_fit(
            power_law,
            age_disp_filtered[age_col],
            age_disp_filtered[disp_col],
            sigma=age_disp_filtered[err_col],
            absolute_sigma=True,
            p0=[1, 0.5]
        )
        sigma_0_fit, beta_fit = popt
        perr = np.sqrt(np.diag(pcov))
        sigma_0_err, beta_err = perr
        age_fit = np.linspace(age_disp_filtered[age_col].min(), age_disp_filtered[age_col].max(), 200)
        ax = axes[i]
        if plot_loglog:
            ax.errorbar(
                age_disp_filtered[age_col], age_disp_filtered[disp_col],
                yerr=age_disp_filtered[err_col],
                fmt='o', color=scatter_color, ecolor='gray', elinewidth=1.5, capsize=3,
                label=f"Mean {labels[comp]} per Age bin"
            )
            ax.plot(
                age_fit, power_law(age_fit, *popt),
                color=fit_color, linestyle=':', linewidth=2,
                label=rf"Fit: $\sigma = {sigma_0_fit:.2f} \times \tau^{{{beta_fit:.2f} \pm {beta_err:.2f}}}$"
            )
            ax.set_xscale('log')
            ax.set_yscale('log')
            # Only set x-axis label for bottom row (indices 2 and 3)
            if i in [2, 3]:
                ax.set_xlabel("Age [Gyr] (log)", fontsize=20)
            else:
                ax.set_xlabel("", fontsize=20)
            ax.set_ylabel(labels[comp] + " (log)", fontsize=20)
            ax.set_title(labels[comp] + " vs Age (log-log)", fontsize=20)
        else:
            ax.errorbar(
                age_disp_filtered[age_col], age_disp_filtered[disp_col],
                yerr=age_disp_filtered[err_col],
                fmt='o', color=scatter_color, ecolor='gray', elinewidth=1.5, capsize=3,
                label=f"Mean {labels[comp]} per Age bin"
            )
            ax.plot(
                age_fit, power_law(age_fit, *popt),
                color=fit_color, linestyle=':', linewidth=2,
                label=rf"Fit: $\sigma = {sigma_0_fit:.2f} \times \tau^{{{beta_fit:.2f} \pm {beta_err:.2f}}}$"
            )
            # Only set x-axis label for bottom row (indices 2 and 3)
            if i in [2, 3]:
                ax.set_xlabel("Age [Gyr]", fontsize=20)
            else:
                ax.set_xlabel("", fontsize=20)
            ax.set_ylabel(labels[comp], fontsize=20)
            ax.set_title(labels[comp] + " vs Age", fontsize=16)
        ax.grid(True, which='both', linestyle='--', alpha=0.3)
        ax.legend(fontsize=14, frameon=True)
    plt.tight_layout()
    plt.savefig(filepath + f"/plots/dispersion/velocity_dispersion_vs_age_allcomponents_withMScuts_pluserrors{'_loglog' if plot_loglog else ''}.png", dpi=300)
    plt.show()

def plot_dispersion_vs_mh(df, filepath, plot_loglog=False):
    mh_col = "param_center"
    disp_col = "sigma_total"
    err_col = "delta_sigma_total"
    # Group by metallicity and compute mean velocity dispersion and error for each metallicity bin
    mh_disp = df.groupby(mh_col)[[disp_col, err_col]].mean().reset_index()
    # Manually specify the bin indices (row numbers) to exclude
    bins_to_exclude = list(range(0, 8)) + list(range(72, 74)) + [60, 81, 88, 96]
    mask = ~mh_disp.index.isin(bins_to_exclude)
    mh_disp_filtered = mh_disp[mask].reset_index(drop=True)

    # Assign disk/halo categories
    def disk_category(mh):
        if mh > -0.5:
            return "Thin Disk"
        elif -1.0 < mh <= -0.5:
            return "Thick Disk"
        elif -1.5 < mh <= -1.0:
            return "Halo"
        else:
            return "Other"
    mh_disp_filtered["disk"] = mh_disp_filtered[mh_col].apply(disk_category)
    disk_colors = {
        "Thin Disk": "#1f77b4",
        "Thick Disk": "#ff7f0e",
        "Halo": "#2ca02c",
        "Other": "#7f7f7f"
    }

    # Option to plot both quadratic and normal (power-law) fit
    plot_normal_fit = True  # Set to True to plot both fits

    # Fit quadratic in log space
    def log_quadratic(mh, a, b, c):
        return a + b * mh + c * mh**2

    popt_quad, pcov_quad = curve_fit(
        log_quadratic,
        mh_disp_filtered[mh_col].values,
        np.log(mh_disp_filtered[disp_col].values),
        sigma=mh_disp_filtered[err_col].values / mh_disp_filtered[disp_col].values,
        absolute_sigma=True,
        p0=[np.log(40), -0.5, 0.1]
    )

    # Fit normal power-law: sigma = sigma_0 * 10^(beta * [M/H])
    def power_law_mh(mh, sigma_0, beta):
        return sigma_0 * np.power(10, beta * mh)

    if plot_normal_fit:
        popt_norm, pcov_norm = curve_fit(
            power_law_mh,
            mh_disp_filtered[mh_col].values,
            mh_disp_filtered[disp_col].values,
            sigma=mh_disp_filtered[err_col].values,
            absolute_sigma=True,
            p0=[40, -0.5]
        )
        sigma_0, beta = popt_norm
        sigma_0_err, beta_err = np.sqrt(np.diag(pcov_norm))

    plt.figure(figsize=(7, 5))
    fit_color_quad = "#d62728"
    fit_color_norm = "#000000"
    mh_fit = np.linspace(mh_disp_filtered[mh_col].min(), mh_disp_filtered[mh_col].max(), 200)
    log_sigma_fit_quad = log_quadratic(mh_fit, *popt_quad)
    sigma_fit_quad = np.exp(log_sigma_fit_quad)
    if plot_normal_fit:
        sigma_fit_norm = power_law_mh(mh_fit, *popt_norm)

    # Plot points by disk category, but do NOT add labels (so legend only shows fits)
    for disk, color in disk_colors.items():
        mask = mh_disp_filtered["disk"] == disk
        if mask.any():
            plt.errorbar(
                mh_disp_filtered.loc[mask, mh_col],
                mh_disp_filtered.loc[mask, disp_col],
                yerr=mh_disp_filtered.loc[mask, err_col],
                fmt='o', color=color, ecolor='gray', elinewidth=1.5, capsize=3,
                label=None  # <-- No label for points
            )

    # Format coefficients with sign
    def fmt_coeff(val):
        return f"{val:.2f}" if val >= 0 else f"- {abs(val):.2f}"

    label_quad = (
        r"Quad Fit: $\log(\sigma) = "
        f"{fmt_coeff(popt_quad[0])} "
        f"{'+' if popt_quad[1] >= 0 else '-'} {abs(popt_quad[1]):.2f} \cdot [\mathrm{{M}}/\mathrm{{H}}] "
        f"{'+' if popt_quad[2] >= 0 else '-'} {abs(popt_quad[2]):.2f} \cdot [\mathrm{{M}}/\mathrm{{H}}]^2$"
    )

    plt.plot(
        mh_fit,
        sigma_fit_quad,
        color=fit_color_quad, linestyle=':', linewidth=2,
        label=label_quad
    )

    if plot_normal_fit:
        label_norm = (
            r"Norm Fit: $\sigma = "
            f"{sigma_0:.2f} \pm {sigma_0_err:.2f} \cdot 10^{{{beta:.2f} \pm {beta_err:.2f} \cdot [\mathrm{{M}}/\mathrm{{H}}]}}$"
        )
        plt.plot(
            mh_fit,
            sigma_fit_norm,
            color=fit_color_norm, linestyle='--', linewidth=2,
            label=label_norm
        )

    if plot_loglog:
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("[M/H] (log)", fontsize=13)
        plt.ylabel(r"$\sigma_\mathrm{total}$ [km/s] (log)", fontsize=13)
        #plt.title(r"Velocity Dispersion vs [M/H] (log-log)", fontsize=15)
    else:
        plt.xlabel("[M/H]", fontsize=13)
        plt.ylabel(r"$\sigma_\mathrm{total}$ [km/s]", fontsize=13)
        #plt.title(r"Velocity Dispersion vs [M/H]", fontsize=15)
    plt.grid(True, which='both', linestyle='--', alpha=0.3)
    plt.legend(fontsize=12, frameon=True)
    plt.tight_layout()
    plt.savefig(filepath + f"/plots/dispersion/velocity_dispersion_vs_mh_bothfits_agecut_maxage{'_loglog' if plot_loglog else ''}.png", dpi=300)
    plt.show()

def plot_powerlaw_trends_vs_age(filepath):
    # Define CSV files and labels
    csv_files = [
        ("combined_sigma_age_binned100.csv", "BaSTI, all stars"),
        ("combined_sigma_age_binned_agecut100.csv", "BaSTI, stars with good ages"),
        ("combined_sigma_age_binned_agecut_verrorcut100.csv", 
         "BaSTI, good ages, $V_\\mathrm{tan} > 100$ km/s, $\\mathrm{error}_{U,V,W} < 2$ km/s"),
        ("combined_sigma_age_binned_age_verror_parallax_MScuts100.csv", 
         "BaSTI, good ages, $V_\\mathrm{tan} > 100$ km/s, $\\mathrm{error}_{U,V,W} < 2$ km/s,\n"
         "$\\sigma_\\mathrm{\\varpi}/\\varpi > 100$, "
         "$T_\\mathrm{eff} > 5500$ K, $L > 0.1\,L_\\odot$")
    ]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    age_col = "param_center"
    disp_col = "sigma_total"
    err_col = "error_sigma_total"
    bins_to_exclude = list(range(1, 8)) + list(range(72, 74)) + [60, 81, 88, 96]
    #bins_to_exclude = None
    def power_law(age, sigma_0, beta):
        return sigma_0 * np.power(age, beta)

    plt.figure(figsize=(8, 6))
    for (csv_name, label), color in zip(csv_files, colors):
        df = pd.read_csv(filepath + f"/outputdata/velocitydispersion/{csv_name}")
        age_disp = df.groupby(age_col)[[disp_col, err_col]].mean().reset_index()
        mask = ~age_disp.index.isin(bins_to_exclude)
        age_disp_filtered = age_disp[mask].reset_index(drop=True)
        popt, pcov = curve_fit(
            power_law,
            age_disp_filtered[age_col],
            age_disp_filtered[disp_col],
            sigma=age_disp_filtered[err_col],
            absolute_sigma=True,
            p0=[1, 0.5]
        )
        age_fit = np.linspace(age_disp_filtered[age_col].min(), age_disp_filtered[age_col].max(), 200)
        plt.plot(
            age_fit, power_law(age_fit, *popt),
            color=color, linestyle='-', linewidth=2,
            label=label  # Only show the cut label, not the formula
        )
    plt.xlabel("Age [Gyr]", fontsize=13)
    plt.ylabel(r"$\sigma_\mathrm{total}$ [km/s]", fontsize=13)
    plt.xlim(0, 14)
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend(fontsize=11, frameon=False)
    plt.tight_layout()
    plt.savefig(filepath + "/plots/dispersion/velocity_dispersion_powerlaw_trends_vs_age.png", dpi=300)
    plt.show()

def main():
    df = load_data(csv_file2, filter_zeros=filter_zeros, components=components)

    #plot_all_components(df, components, labels, cmap, filepath, parameter_names)
    
    #plot_single_component(df, single_component, labels, cmap, filepath)
    #plot_dispersion_vs_age(df, filepath, plot_loglog=False)
    #plot_dispersion_vs_age(df, filepath, plot_loglog=False)
    plot_dispersion_vs_age_all_components(df, filepath, plot_loglog=False)
    #plot_dispersion_vs_mh(df, filepath, plot_loglog=False)
    #plot_powerlaw_trends_vs_age(filepath)
    
if __name__ == "__main__":
    main()
