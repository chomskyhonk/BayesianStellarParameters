import matplotlib.pyplot as plt
import numpy as np
import os

# File containing the isochrone data
file_path = "/Users/sarahharrison/projects/vsc/iniplots/t449_ACS_z001_y0259/isochrone_datasets/isochrones/isochrone_age_100_Z_z01.txt"

# Initialize lists to store isochrone data
logTe_lines = []
logLLo_lines = []

# Read data from the file
with open(file_path, "r") as file:
    for line in file:
        values = line.strip().split()
        if len(values) == 2:
            logTe, logLLo = map(float, values)
            logTe_lines.append(logTe)
            logLLo_lines.append(logLLo)

# Star data
starLogTe = [3.75]
starLogLLo = [0.4]

# Number of points per isochrone (i.e., break the line every 2000 points)
chunk_size = 2000

# Create a figure and axis for the plot
fig, ax = plt.subplots()

# Loop through the isochrone data in chunks and plot each chunk separately
for i in range(1000, len(logTe_lines), chunk_size):
    # Get the current chunk of data
    chunk_logTe = logTe_lines[i:i + chunk_size]
    chunk_logLLo = logLLo_lines[i:i + chunk_size]
    
    # Plot the chunk as a line
    ax.plot(chunk_logTe, chunk_logLLo, color='black', linewidth=0.3, zorder=1)
    
    # Optional: If this is not the last chunk, connect it to the next chunk with a line
    #if i + chunk_size < len(logTe_lines):
    #    ax.plot([chunk_logTe[-1], logTe_lines[i + chunk_size]], 
    #            [chunk_logLLo[-1], logLLo_lines[i + chunk_size]], 
    #            color='red', linewidth=0.3, zorder=1)

# Plot scatter points on top
ax.scatter(starLogTe, starLogLLo, s=50, color='red', zorder=2)

# Add labels, limits, and title
ax.set_xlabel("logTe")
ax.set_xlim(4.3, 3.4)  # Reverse the x-axis limits
ax.set_ylim(-2, 5)
ax.set_ylabel("Log(L/Lo)")
ax.set_title("Luminosity vs T_eff")
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "isochronegrid.png")
plt.savefig(output_path)

# Show the plot
plt.show()