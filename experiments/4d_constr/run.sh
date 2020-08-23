#!/bin/bash

declare -a threads=("1" "2" "4" "6")
declare -a delays=("0.5" "1")
declare -a modes=("multi" "asynch")
executable='../bin/solve_constrained_sample'

dim="4"

for mode in "${modes[@]}"
do
  mode_str=${mode//[.]/}
  for delay in "${delays[@]}"
  do
    delay_str=${delay//[.]/}
    for t_num in "${threads[@]}"
    do
      command="$executable -n 20 -p $t_num -m $mode -d $dim -e 0.02 -r 4.7 -s --delay $delay -f gkls_"$dim"d_m_"$mode_str"_p_"$t_num"_d_$delay_str.csv"
      echo $command
      $command
    done
  done
done
