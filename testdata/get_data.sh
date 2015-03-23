#! /bin/bash
RESTDIR="/scratch/braindata/eglerean/rest_vs_task_HO/rest/"
TASKDIR="/scratch/braindata/eglerean/rest_vs_task_HO/task/"

for i in `seq 1 13`;
    do
        cp $RESTDIR$i"/net.mat" "rest_HO_"$i".mat"
        cp $TASKDIR$i"/net.mat" "movie_HO_"$i".mat"
    done
