#!/bin/bash

#echo 'debug version'
#time ../../debug/bin/parallel_track_example lat.bmad


export OMP_NUM_THREADS='10'
echo '10 threads'
time ../../production/bin/parallel_track_example lat.bmad

export OMP_NUM_THREADS=1
echo '1 thread'
time ../../production/bin/parallel_track_example lat.bmad
