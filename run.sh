#!/bin/bash
mkdir build
cd build 
cmake ../
make -j16
./MP
python3 ../animate.py
python3 ../threads.py
