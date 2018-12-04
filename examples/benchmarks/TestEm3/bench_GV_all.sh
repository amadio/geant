#!/bin/bash
# F=field
# P=physics
# G=geometry
# M=MSC
#___________________F_P_G_M
./TestEm3_GV.run "$1" 0 0 0 0
./TestEm3_GV.run "$1" 1 0 0 0
./TestEm3_GV.run "$1" 2 0 0 0
./TestEm3_GV.run "$1" 0 1 0 0
./TestEm3_GV.run "$1" 0 2 0 0
./TestEm3_GV.run "$1" 0 0 1 0
./TestEm3_GV.run "$1" 0 0 2 0
./TestEm3_GV.run "$1" 0 0 0 1
./TestEm3_GV.run "$1" 0 0 0 2
