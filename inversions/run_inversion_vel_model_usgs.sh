#!/bin/bash
set -euo pipefail

#lon=124.5
lon=124.75
myfolder=usgs_vel

# Prompt for user input
echo "Do you want to rerun the gps (y/n)"
read -r rerun_gps

echo "Do you want to rerun the gps-teleseismic inversion (y/n)"
read -r rerun_gps_tele

echo "Do you want to rerun the gps-teleseismic-strong motion inversion (y/n)"
read -r rerun_gps_tele_strong

suffix=gps
if [[ "$rerun_gps" == "y" ]]; then
  #we first run a auto inversion
  ffm model run $(pwd) auto_model -g data/cmtsolution -t gps -d data/Static_Data/

  cp 20241205184419/ffm.0/NP1/ ${myfolder}_$suffix -r
  rm -r 20241205184419

  #python scripts/prepare_velocity_model_canvas.py $lon
  #python scripts/prepare_velocity_model.py $lon
  cp input_data/segments_data.json ${myfolder}_$suffix
  cp input_data/annealing_prop.json ${myfolder}_$suffix
  #cp data/vel_model.txt ${myfolder}_$suffix

  cd ${myfolder}_${suffix}
  #ffm model run $(pwd) manual_model_add_data
  #ffm manage velmodel-to-json $(pwd) vel_model.txt
  ffm manage update-inputs $(pwd) -p -m -a
  ffm model run $(pwd) manual_model_add_data
  cp Solution.txt plots
  #cp modelling_summary.txt plots
  cd ..
fi

if [[ "$rerun_gps_tele" == "y" ]]; then
  cp -r ${myfolder}_${suffix} ${myfolder}_gps_tele
  suffix=gps_tele
  cd ${myfolder}_${suffix}

  echo "starting inv with teleseismics"
  ffm model run $(pwd) manual_model_add_data -t body -t surf -d ../data/Teleseismic_Data/
  ffm manage modify-dicts $(pwd) downweight surf -sc "DWPF:BHZ,SH"
  ffm manage modify-dicts $(pwd) downweight surf -sc "MIDW:BHT"
  ffm manage modify-dicts $(pwd) downweight surf -sc "HRV:BHT"
  ffm manage update-inputs $(pwd) -t surf

  cp Solution.txt plots
  cp modelling_summary.txt plots
  mv plots plots_gps_tele

  ffm process shift-match $(pwd) body -o auto
  ffm process shift-match $(pwd) surf -o auto
  ffm process remove-baseline $(pwd)
  ffm manage update-inputs $(pwd) -t body -t surf
  echo "starting inv with teleseismics after shift-match"
  ffm model run $(pwd) manual_model_add_data
  cp Solution.txt plots
  cp modelling_summary.txt plots
  mv plots plots_gps_tele_shift_match
  cd ..
fi

if [[ "$rerun_gps_tele_strong" == "y" ]]; then
  suffix=gps_tele
  cp -r ${myfolder}_$suffix ${myfolder}_gps_tele_sm
  suffix=gps_tele_sm
  cd ${myfolder}_${suffix}
  ffm model run $(pwd) manual_model_add_data -t strong -d ../data/StrongMotion_Data/
  cp Solution.txt plots
  mv plots plots_gps_tele_sm
  echo "done running inversion with strong motion"

  ../../submodules/seismic-waveform-factory/scripts/modify_wasp_strong_motion_waves.py ../input_data/waveforms_config.ini
  ffm model run $(pwd) manual_model_add_data
  cp Solution.txt plots
  cp modelling_summary.txt plots
  mv plots plots_gps_tele_sm_shorter
  cd ..
  echo "done running inversion with strong motion (shorter)"

  suffix=gps_tele_sm
  cp -r ${myfolder}_${suffix} ${myfolder}_gps_tl_sm2
  suffix=gps_tl_sm2
  cd ${myfolder}_${suffix}
  ffm process shift-match $(pwd) strong -o auto
  ffm model run $(pwd) manual_model_add_data
  cp Solution.txt plots
  cp modelling_summary.txt plots
  cd ..
  echo "done running inversion with strong motion (shorter + shift-match)"
fi
