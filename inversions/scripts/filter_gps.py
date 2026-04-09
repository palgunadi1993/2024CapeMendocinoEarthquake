#!/usr/bin/env python3
import pandas as pd
import numpy as np

# Load the file (adjust the filename as needed)
df = pd.read_csv("data/nc75095651_web.txt", sep='\s+', skiprows=[1])
print(df.keys())

filtered_df = df[(df["Lon"] < -123.2) & (df["Lat"] > 39.75) & (df["Lat"] < 41.5)]

# Define the output file
output_file = "data/Static_Data/gnss_data"

# Save to file with the same formatting
with open(output_file, "w") as f:
    # Write the header
    f.write(
        "Sta    Lon       Lat       de(m)     dn(m)     du(m)    sde(m)   sdn(m)   sdu(m)\n"
    )
    f.write("=" * 80 + "\n")  # Add the separator line

    # Save data without the index, using fixed-width formatting
    filtered_df.to_string(f, index=False, header=False, justify="left")

print(f"Filtered data saved to {output_file}")
