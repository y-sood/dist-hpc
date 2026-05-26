import matplotlib.pyplot as plt

# Read data from results.txt
with open('results.txt', 'r') as f:
    data = [float(line.strip()) for line in f if line.strip()]

# Plot histogram
plt.hist(data, bins=20, edgecolor='black')
plt.xlabel('Number of Susceptible Humans')
plt.ylabel('Frequency')
plt.title('Histogram of Susceptible Humans after Monte Carlo Simulation')
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('histogram.png', dpi=300)