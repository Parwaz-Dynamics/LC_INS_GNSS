import pandas as pd
from collections import Counter

# File path as provided
file_path = r'C:\Users\HP\Desktop\Abdullah Wasim\Datasets\GPS.csv'

total = sum(1 for _ in open(file_path)) - 1
half = total

chunks = []
count = 0
for chunk in pd.read_csv(file_path, chunksize=100000):
    chunks.append(chunk)
    count += len(chunk)
    if count >= half:
        break
df = pd.concat(chunks)[:half]

# Extract the first column (gps_obssec(s))
# Assuming the first column is index 0, and it may have a header
times = df.iloc[:, 0].values

# If the CSV has headers, and the column name is 'gps_obssec(s)', use:
# times = df['gps_obssec(s)'].values

# Compute differences dt = t(k) - t(k-1)
dts = [times[i] - times[i-1] for i in range(1, len(times))]

# To handle floating-point precision, round dts to a certain decimal places (e.g., 6)
rounded_dts = [round(dt, 6) for dt in dts]

# Count the frequency of each dt value
dt_counts = Counter(rounded_dts)

# Sort by dt value and print the results
print("Analysis of dt values and their frequencies:")
for dt, count in sorted(dt_counts.items()):
    print(f"dt: {dt}, count: {count}")

# Find indices where 0.374, 0.636, and 5.47 occur in rounded_dts
targets = [37.005, 0.048, 0.165]
for target in targets:
    indices = [i for i, val in enumerate(rounded_dts) if val == target]
    # print(f"Indices where dt == {target}: {times[indices]}")