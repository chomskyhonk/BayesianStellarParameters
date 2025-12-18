import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

data = pd.read_csv(
    "/path/to/your/data",
    on_bad_lines='skip'
)

age_bins = np.linspace(0.1, 14, 20)
mh_bins = np.linspace(-1.5, 0.42, 20)

data['age_bin'] = pd.cut(data['age_flame'], bins=age_bins, labels=False, include_lowest=True)
data['mh_bin'] = pd.cut(data['mh_gspspec'], bins=mh_bins, labels=False, include_lowest=True)
#replace age_flame and mh_gspsec with your derived Bayesian ages/metallicity columns
results = []

for age_bin in range(len(age_bins)-1):
    for mh_bin in range(len(mh_bins)-1):
        subset = data[(data['age_bin'] == age_bin) & (data['mh_bin'] == mh_bin)]
        if len(subset) > 0:
            u_disp = np.std(subset['U'])
            v_disp = np.std(subset['V'])
            w_disp = np.std(subset['W'])
            total_disp = np.sqrt(u_disp**2 + v_disp**2 + w_disp**2)
            results.append({
                'age_bin': age_bin,
                'mh_bin': mh_bin,
                'u_disp': u_disp,
                'v_disp': v_disp,
                'w_disp': w_disp,
                'total_disp': total_disp,
                'n_stars': len(subset)
            })

disp_df = pd.DataFrame(results)

#discard bins where total_disp is 0 or NaN, or n_stars < 2
disp_df = disp_df[(disp_df['total_disp'] > 0) & disp_df['total_disp'].notna() & (disp_df['n_stars'] >= 10)]

data['Ulsr'] = data['U'] + 11.10
data['Vlsr'] = data['V'] + 12.24 - 232
data['Wlsr'] = data['W'] + 7.25
data['V_tot'] = np.sqrt(data['Ulsr']**2 + data['Vlsr']**2 + data['Wlsr']**2)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

#compute average V_tot for each bin
avg_vtot = []
for _, row in disp_df.iterrows():
    subset = data[(data['age_bin'] == row['age_bin']) & (data['mh_bin'] == row['mh_bin'])]
    avg_vtot.append(subset['V_tot'].mean() if len(subset) > 0 else np.nan)
disp_df['avg_vtot'] = avg_vtot

#bin centers for plotting
age_centers = 0.5 * (age_bins[:-1] + age_bins[1:])
mh_centers = 0.5 * (mh_bins[:-1] + mh_bins[1:])

#map bin indices to bin centers
disp_df['age_center'] = disp_df['age_bin'].apply(lambda x: age_centers[int(x)])
disp_df['mh_center'] = disp_df['mh_bin'].apply(lambda x: mh_centers[int(x)])
sc = ax.scatter(
    disp_df['age_center'],
    disp_df['mh_center'],
    disp_df['total_disp'],
    c=disp_df['avg_vtot'],
    cmap='viridis',
    s=60
)

ax.set_xlabel('Age [Gyr]', fontsize=18)
ax.set_ylabel('[M/H]', fontsize=18)
ax.set_zlabel(r'$\sigma_{\mathrm{tot}}$ [km/s]', fontsize=18)

cb = fig.colorbar(sc, ax=ax, pad=0.1, shrink=0.9, orientation='vertical')
#cb.set_label(r'Average $V_{\mathrm{tot}}$ [km/s]', labelpad=10)
cb.ax.set_xlabel(r'Average $V_{\mathrm{tot}}$ [km/s]', labelpad=10, fontsize=16)
cb.ax.xaxis.set_label_position('bottom')

plt.tight_layout()
plt.savefig("/path/to/your/savefigfile")
plt.show()

