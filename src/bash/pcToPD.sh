#! /bin/bash

FILES=/home/elove/Desktop/actin_07102019/sampled_point_clouds/*

for file  in $FILES 
do
  /home/elove/Desktop/ripser/ripser --threshold 10 --dim 1 --format point-cloud $file > "/home/elove/Desktop/actin_07102019/sampled_pds/"${file##/*/}
done
