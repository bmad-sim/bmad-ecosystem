import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def read_turn_file(filename):
    with open(filename, 'r') as f:
        # Skip header lines until we find the line starting with ##
        for line in f:
            if line.startswith('##'):
                # Remove '##', strip whitespace, and split by '|'
                header_parts = line.strip('#').strip().split('|')
                # Split by whitespace to get column names
                header = ' '.join(header_parts).split()
                break
        
        # Read the data into a DataFrame
        df = pd.read_csv(f, sep='\s+', names=header)
    
    # Extract turn number from filename
    turn = int(filename.split('turn')[1])
    df['Turn'] = turn
    
    return df

# Get all turn files in the current directory
turn_files = sorted(glob.glob('turn*'))

# Read all turn files
all_data = pd.concat([read_turn_file(file) for file in turn_files], ignore_index=True)

# Create the scatter plot
plt.figure(figsize=(12, 8))

# Get unique particle IDs
particle_ids = all_data['Ix'].unique()

# Plot trajectory for each particle
for particle_id in particle_ids:
    particle_data = all_data[all_data['Ix'] == particle_id]
    plt.scatter(particle_data['z'], particle_data['pz'], s=10, label=f'Particle {particle_id}')

plt.xlabel('z')
plt.ylabel('pz')
plt.title('Particle Trajectories in z-pz Plane')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()