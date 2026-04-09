#!/bin/bash
set -euo pipefail

# Prompt for user input
echo "Do you want to download CMT data? (y/n)"
read -r download_cmt

echo "Do you want to download GPS data? (y/n)"
read -r download_gps

echo "Do you want to download Teleseismic data? (y/n)"
read -r download_teleseismic

echo "Do you want to download Strong Motion data? (y/n)"
read -r download_strongmotion

echo "Do you want to prepare the velocity model? (y/n)"
read -r prepare_velocity_model

# Conditional execution based on user input

# Download CMT if needed
if [[ "$download_cmt" == "y" || "$download_cmt" == "Y" ]]; then
  echo "Downloading CMT data from GlobalCMT..."
  mkdir -p data
  scripts/download_cmtsolution.py --date 2024-12-05 --lat 40.35 --lon -125.0 --min-mag 6.5 -o data/cmtsolution
  scripts/modify_depth_in_cmtsolution.py
else
  echo "Skipping CMT download."
fi

# Download GPS if needed
if [[ "$download_gps" == "y" || "$download_gps" == "Y" ]]; then
  echo "Downloading GPS data..."
  wget https://geodesy.unr.edu/news_items/20241211/nc75095651_web.txt --directory-prefix data
  mkdir -p data/Static_Data
  scripts/filter_gps.py
else
  echo "Skipping GPS download."
fi

# Download Teleseismic data if needed
if [[ "$download_teleseismic" == "y" || "$download_teleseismic" == "Y" ]]; then
  echo "Downloading Teleseismic data..."
  mkdir -p data/Teleseismic_Data
  ffm manage acquire data/Teleseismic_Data data/cmtsolution -t body
else
  echo "Skipping Teleseismic download."
fi

# Download Strong Motion data if needed
if [[ "$download_strongmotion" == "y" || "$download_strongmotion" == "Y" ]]; then
  echo "for missing modules, consider"
  echo "pip install -r ../submodules/seismic-waveform-factory/requirements.txt"
  echo "Downloading Strong Motion data..."
  mkdir -p data/StrongMotion_Data
  ../submodules/seismic-waveform-factory/scripts/select_stations.py input_data/waveforms_config.ini 30 10 --channel "HN*"
  ../submodules/seismic-waveform-factory/scripts/copy_selected_sac_data_to_folder.py input_data/waveforms_config.ini --output_folder data/StrongMotion_Data/
else
  echo "Skipping Strong Motion download."
fi

# Download Strong Motion data if needed
if [[ "$prepare_velocity_model" == "y" || "$prepare_velocity_model" == "Y" ]]; then

  FILE="data/casc1.6-velmdl.r1.0-n4.nc"
  echo "Downloading casc1.6..."
  echo "this will download the 5km sampled model. For the full resolution model, go to https://www.sciencebase.gov/catalog/item/59f1e68be4b0220bbd9dd4b4"

  wget -O "$FILE" "https://cvm.cascadiaquakes.org/data/download-netcdf-s3/?filename=casc1.6-velmdl.r1.0-n4.json"

  echo "Velocity model downloaded. Run prepare_velocity_model.py with a longitude argument during inversion (e.g. scripts/prepare_velocity_model.py 124.75)"
else
  echo "Skipping preparing velocity model"
fi
