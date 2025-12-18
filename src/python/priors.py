import numpy as np
import matplotlib.pyplot as plt


def metallicity_prior(mh, prior_type="gaussian"):
    """
    Calculate the metallicity prior based on the given M/H value.
    Options:
    - "flat": Flat prior (constant value of 1.0).
    - "gaussian": Gaussian prior peaked at -0.1 dex with uncertainty of 0.3 dex.
    """
    if prior_type == "flat":
        return 1.0  # Flat prior
    elif prior_type == "gaussian":
        mean = -0.250103
        sigma = 0.423109
        return np.exp(-0.5 * ((mh - mean) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    else:
        raise ValueError(f"Unknown prior type: {prior_type}")


def age_prior(age, mh):
    """
    Calculate the age prior based on the given age and M/H value.
    """
    if age > 14.0:
        return 0.0
    elif 11.0 <= age <= 14.0:
        return 1.0
    else:
        # Determine sigmat based on M/H
        if mh < -0.9:
            sigmat = 1.5
        elif -0.9 <= mh <= -0.5:
            sigmat = 1.5 + 7.5 * (0.9 + mh) / 0.4
        else:
            sigmat = 9.0
        return np.exp((age - 11.0) / sigmat)

def total_prior(age, mh):
    """
    Calculate the total prior as the product of the age prior and the metallicity prior.
    """
    return age_prior(age, mh) * metallicity_prior(mh)

def save_prior_values_to_text(prior_values, filename):
    """
    Save the prior values to a text file.

    Args:
        prior_values (np.ndarray): 2D array of prior values.
        filename (str): Path to the output text file.
    """
    with open(filename, 'w') as f:
        for row in prior_values:
            f.write(" ".join(map(str, row)) + "\n")
    print(f"Prior values saved to {filename}")

def plot_combined_prior(min_age, max_age, min_mh, max_mh, num_age_bins, num_mh_bins, savelocation="/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/plots/combined_prior.png", logscale=True):    
    """
    Plot the combined prior function for the given ranges and bins.
    """
    # Create age and metallicity bins
    age_bins = np.linspace(min_age, max_age, num_age_bins)
    mh_bins = np.linspace(min_mh, max_mh, num_mh_bins)

    # Initialize a 2D array to store the prior values
    prior_values = np.zeros((num_age_bins, num_mh_bins))

    # Calculate the prior for each bin
    for i, age in enumerate(age_bins):
        for j, mh in enumerate(mh_bins):
            prior_values[i, j] = total_prior(age, mh)
            
    # Normalize the prior values so they sum to 1
    total_sum = np.sum(prior_values)
    if total_sum > 0:
        prior_values /= total_sum

    # Optional: log-scale transform for visualization
    if logscale:
        with np.errstate(divide='ignore'):
            prior_values = np.log10(prior_values)
            # Set -inf (log(0)) to minimum finite log value for color scaling
            finite_mask = np.isfinite(prior_values)
            min_val = np.min(prior_values[finite_mask]) if np.any(finite_mask) else 0
            prior_values[~finite_mask] = min_val

    # Save to file
    output_file = "/path/to/your/files/outputdata/prior_values.txt"
    save_prior_values_to_text(prior_values, output_file)

    # Plot the combined prior as a heatmap
    plt.figure(figsize=(8, 6))
    extent = [min_mh, max_mh, min_age, max_age]
    im = plt.imshow(prior_values, extent=extent, origin='lower', aspect='auto', cmap='viridis')
    cbar = plt.colorbar(im, orientation='vertical', pad=0.04)
    cbar.ax.tick_params(labelsize=13)
    # Remove the default label
    cbar.set_label('')
    # Add a label below the colorbar, aligned with the x-axis
    cbar.ax.set_xlabel('P(Age, [M/H])', labelpad=15, fontsize=13)
    cbar.ax.tick_params(labelsize=13)
    plt.xlabel('Metallicity [M/H]', fontsize=13)
    plt.ylabel('Age [Gyr]', fontsize=13)
    #plt.title('Combined Prior Function', fontsize=13)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.savefig(savelocation)
    plt.show()
        


# Example usage
if __name__ == "__main__":
    # Define the ranges and parameters
    min_age = 0.1  # Minimum age (Gyr)
    max_age = 14.0  # Maximum age (Gyr)
    min_mh = -1.5  # Minimum metallicity [M/H]
    max_mh = 0.42  # Maximum metallicity [M/H]
    num_age_bins = 50  # Number of age bins
    num_mh_bins = 48  # Number of metallicity bins

    # Plot the combined prior
    # Plot the combined prior with linear scale
    # plot_combined_prior(min_age, max_age, min_mh, max_mh, num_age_bins, num_mh_bins, logscale=False)
    # Plot the combined prior with log scale
    plot_combined_prior(min_age, max_age, min_mh, max_mh, num_age_bins, num_mh_bins, logscale=False)
