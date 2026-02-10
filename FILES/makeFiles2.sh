mkdir ./FILES
cd ./FILES
for i in $(seq 1 100); do touch file${i}.txt; done
