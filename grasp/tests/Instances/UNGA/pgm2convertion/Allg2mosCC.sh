#!/bin/bash

./g2mosCC ../UNGA_Weighted_SignedGraphs/Section01.g
./g2mosCC ../UNGA_Weighted_SignedGraphs/Section02.g
./g2mosCC ../UNGA_Weighted_SignedGraphs/Section03.g
./g2mosCC ../UNGA_Weighted_SignedGraphs/Section04.g
./g2mosCC ../UNGA_Weighted_SignedGraphs/Section05.g
./g2mosCC ../UNGA_Weighted_SignedGraphs/Section06.g
./g2mosCC ../UNGA_Weighted_SignedGraphs/Section07.g
./g2mosCC ../UNGA_Weighted_SignedGraphs/Section08.g
./g2mosCC ../UNGA_Weighted_SignedGraphs/Section09.g

for i in {10..63}
do
   ./g2mosCC ../UNGA_Weighted_SignedGraphs/Section$i.g
done
