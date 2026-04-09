#!/usr/bin/env python3
from netCDF4 import Dataset
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description="prepare a velocity model from ffm, extracting from casc1.6 model"
)
parser.add_argument("lon", type=float, help="reference longitude (in absolute value)")
args = parser.parse_args()


# File path
file_path = "data/casc1.6-velmdl.r1.0-n4.nc"

# Open the NetCDF file
with Dataset(file_path, mode="r") as nc_file:
    # Read variables
    lon = nc_file.variables["longitude"][:]
    lat = nc_file.variables["latitude"][:]
    depth = nc_file.variables["depth"][:].astype(float)
    Vp = nc_file.variables["Vp"][:]  # P-wave velocity
    Vs = nc_file.variables["Vs"][:]  # S-wave velocity


extract_single_profile = True

if extract_single_profile:
    target_lon, target_lat = -args.lon, 40.5
    print(target_lon, target_lat)
    # target_lon, target_lat = -124, 40.5
    distances = np.sqrt(np.power(lat - target_lat, 2) + np.power(lon - target_lon, 2))
    closest_index = np.argmin(distances)
    # Convert the flat index to 2D indices (row, col)
    closest_index_2d = np.unravel_index(closest_index, distances.shape)
    i, j = closest_index_2d

    # remove water layer
    lateral_avg_Vp = Vp[:, i, j]
    lateral_avg_Vs = Vs[:, i, j]
    first_non_zero = np.where(lateral_avg_Vs > 0)[0][0]
    lateral_avg_Vp[0:first_non_zero] = lateral_avg_Vp[first_non_zero]
    lateral_avg_Vs[0:first_non_zero] = lateral_avg_Vs[first_non_zero]

else:
    # compute average velocity model over a region
    # Define the area of interest
    # lon_min, lon_max = -125.4, -124
    # lat_min, lat_max = 40.2, 40.5
    # lon_min, lon_max = -125.4, -125.2
    # lat_min, lat_max = 40.2, 40.4

    # lon_min, lon_max = -125.4, -123
    # lat_min, lat_max = 39.5, 41.5
    lon_min, lon_max = -125.4, -124
    lat_min, lat_max = 40.2, 40.5

    # Create a mask for the area of interest
    mask = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)

    # Loop over depth levels and compute the lateral average
    lateral_avg_Vp = []
    lateral_avg_Vs = []

    for d in range(len(depth)):
        # Mask Vp and Vs at the current depth
        Vp_depth = np.ma.masked_array(Vp[d, :, :], mask=~mask)
        Vs_depth = np.ma.masked_array(Vs[d, :, :], mask=~mask)

        # Compute the mean excluding NaN values
        avg_Vp = np.median(Vp_depth[~Vp_depth.mask])
        avg_Vs = np.median(Vs_depth[~Vs_depth.mask])

        lateral_avg_Vp.append(avg_Vp)
        lateral_avg_Vs.append(avg_Vs)

    # Convert to arrays for easy handling
    lateral_avg_Vp = np.array(lateral_avg_Vp)
    lateral_avg_Vs = np.array(lateral_avg_Vs)

# cut at mantle
first_id_mantle = np.where(lateral_avg_Vp > 8000.0)[0][0]
lateral_avg_Vp = lateral_avg_Vp[0:first_id_mantle]
lateral_avg_Vs = lateral_avg_Vs[0:first_id_mantle]
depth = depth[0:first_id_mantle]

# Print results
print("Depth levels (m):", depth)
print("Lateral Average Vp (m/s):", lateral_avg_Vp)
print("Lateral Average Vs (m/s):", lateral_avg_Vs)


plot_velocity_model = False
if plot_velocity_model:
    import matplotlib.pyplot as plt
    import numpy as np

    # Create the plot
    plt.figure(figsize=(8, 6))

    # Plot Vp
    plt.plot(
        lateral_avg_Vp, depth, label="Vp (P-wave velocity)", color="blue", marker="o"
    )
    # Plot Vs
    plt.plot(
        lateral_avg_Vs, depth, label="Vs (S-wave velocity)", color="red", marker="x"
    )

    # Customize the plot
    plt.gca().invert_yaxis()  # Invert y-axis so depth increases downward
    plt.xlabel("Velocity (m/s)", fontsize=12)
    plt.ylabel("Depth (m)", fontsize=12)
    plt.title("Lateral Average Velocity vs Depth", fontsize=14)
    plt.legend()
    plt.grid(True)

    # Show the plot
    plt.tight_layout()
    plt.savefig("average_velocity_model.png")


depth_layer = np.diff(depth) / 1000.0
nlayer = depth_layer.size
lateral_avg_Vp /= 1000.0
lateral_avg_Vs /= 1000.0


def compute_rho(vp_utm):
    # Compute rho from vp
    rho = (
        1.6612 * vp_utm
        - 0.4721 * vp_utm**2
        + 0.0671 * vp_utm**3
        - 0.0043 * vp_utm**4
        + 0.000106 * vp_utm**5
    )
    return rho


rho = compute_rho(lateral_avg_Vp)
# rho = np.linspace(0, nlayer, nlayer+1)


def reduce(rho):
    return [*rho[0:4], *rho[4:10:2], *rho[10::4]]


reduce_nb_layers = False
if not reduce_nb_layers:
    dens, vp, vs, qa, qb, thick = [
        val.data
        for val in [
            rho,
            lateral_avg_Vp,
            lateral_avg_Vs,
            2000 * lateral_avg_Vs,
            1000 * lateral_avg_Vs,
            depth_layer,
        ]
    ]
else:
    depth_layer_less = [
        *depth_layer[0:4],
        *depth_layer[4:10:2] * 2,
        *depth_layer[10::4] * 4,
    ]
    # depth_layer_less = [str(x) for x in depth_layer_less]

    dens = reduce(rho)
    vp = reduce(lateral_avg_Vp)
    vs = reduce(lateral_avg_Vs)
    qa = reduce(2000 * lateral_avg_Vs)
    qb = reduce(1000 * lateral_avg_Vs)
    thick = depth_layer_less

data = {"dens": dens, "p_vel": vp, "s_vel": vs, "thick": thick, "qa": qa, "qb": qb}
print(data)

n = len(thick)
out = f"{n}\n"
for i in range(n):
    out += f"{vp[i]} {vs[i]} {rho[i]} {thick[i]} {qa[i]} {qb[i]}\n"

depth_max_model = np.array(thick).sum()
print(depth_max_model)
depth_first_layer_mantle = 196 - (depth_max_model - 20)

mantle = f"""8.08 4.473 3.3754 {depth_first_layer_mantle} 1200.0 500.0
8.594 4.657 3.4465 36.0 360.0 140.0"""

mantle_axitra = f"""{depth_first_layer_mantle} 8.08 4.473 3.3754 1200.0 500.0
10e3 8.594 4.657 3.4465 360.0 140.0"""

with open("data/vel_model.txt", "w") as fid:
    fid.write(out)
    fid.write(mantle)

# thick[-1] = 10000.0
out = "H P_VEL S_VEL DENS QP QS\n"
for i in range(n):
    out += f"{thick[i]} {vp[i]} {vs[i]} {rho[i]} {qa[i]} {qb[i]}\n"

out += mantle_axitra

fn = "data/vel_model_axitra_fmt.txt"
with open(fn, "w") as fid:
    fid.write(out)
print(f"done writing {fn}")
