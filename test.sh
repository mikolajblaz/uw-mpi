
#!/bin/bash

for i in 1 2 3
do
  time mpirun -np 4 ./collisions-$i -v --hor 2 --ver 2 --gal1 gal1.txt --gal2 gal2.txt --delta 0.1 --total 1.0
done
