#!/bin/bash

# for (\Delta t)^4 fit

# for dt in $(seq 0.01 0.01 0.7)
for dt in $(seq 0.001 0.001 0.009)
do
    sed -i "s/dt = .* # time step/dt = $(echo $dt) # time step/" Ex2_RK4_variables.py
    python3 error.py >> global_ERR
    sleep 0.8 # I don't understand why python would not import the updated file Ex2_RK4_variables.py but a version shortly before instead, if the program is not paused.
done


# for ERR.h5
# for dt in $(seq 0.05 0.05 0.4) 0.5
for dt in $(seq 0.01 0.01 0.7)
do
    sed -i "s/dt = .* # time step/dt = $(echo $dt) # time step/" Ex2_RK4_variables.py
    python3 error.py 1
    sleep 0.8 # I don't understand why python would not import the updated file Ex2_RK4_variables.py but a version shortly before instead, if the program is not paused.
done
