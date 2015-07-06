#!/bin/bash

./mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --branch $2 -N  5 -K 4 < in.0 > out-04-05.txt
./mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --branch $2 -N  6 -K 4 < in.0 > out-04-06.txt
./mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --branch $2 -N  7 -K 4 < in.0 > out-04-07.txt
./mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --branch $2 -N  9 -K 4 < in.0 > out-04-09.txt
./mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --branch $2 -N 10 -K 4 < in.0 > out-04-10.txt
./mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --branch $2 -N 11 -K 4 < in.0 > out-04-11.txt
./mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --branch $2 -N 13 -K 4 < in.0 > out-04-13.txt
./mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --branch $2 -N 14 -K 4 < in.0 > out-04-14.txt
./mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --branch $2 -N 15 -K 4 < in.0 > out-04-15.txt
./mmsconj.exe run -m 1000 -k 3600 --maxsols 10 --prop $1 --branch $2 -N 17 -K 4 < in.0 > out-04-17.txt 
