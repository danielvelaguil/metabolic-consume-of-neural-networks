#!/usr/bin/env python
# -*- coding: utf-8 -*-


#  file avalanches_git.py  Daniel Aguilar Velazquez danielvelaguil@gmail.com
#     This script belongs to the paper "Low metabolic cost of rich-club neural networks at criticality"
# The code print the total number of spikes and synaptic cost of avalanche activity (metabolic_cost_neural_network.cpp) in a 
# hierarchical neural network. An example of the use of this code in ubuntu terminal:
#  ./metabolic_cost_neural_network.x 8.0 .9 .75 100000 7 2.5 1.5 |./avalanches_git.py


import string
import sys, math, os, shutil
import numpy as np

lines=sys.stdin.readlines()
spikes=[]
synaptic_cost=0.0
 
def search(spikes):
 list_n=[]
 for i in spikes:
  if(len(list_n)==0):
   list_n.append(i)
  else:
   flag=0
   for j in list_n:
    if(i==j):
     flag=1
     break;
   if(flag==0):
    list_n.append(i)
 return len(list_n)

for line in lines:
 c1,c2=string.split(line)
 if(float(c1)==-1):
  print search(spikes), len(spikes), synaptic_cost # column 1 prints the size of the avalanche, column 2 prints the total
  spikes=[]                      #number of spikes, column 3 prints the synaptic cost 
  synaptic_cost=0.0
 else:
  if(float(c1)!=0):
   spikes.append(float(c1))
  synaptic_cost+=float(c2)
