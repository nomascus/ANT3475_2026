mkdir ./FILES
cd ./FILES
for i in $(seq 1 50); do touch file${i}.txt; done
