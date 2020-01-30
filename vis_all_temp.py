#!/usr/bin/env python

from athena_read import *
from numpy import *
from matplotlib.pyplot import *
from scipy import stats
from argparse import *
from vis_temp import *

parser = ArgumentParser()

parser.add_argument('--prob')
parser.add_argument('--out')
parser.add_argument('--num')
parser.add_argument('--xmax')
parser.add_argument('--ymax')
parser.add_argument('--qmax')
parser.add_argument('--fnm')
parser.add_argument('-v',action='store_true',default=False)

args = parser.parse_args()

prob = args.prob
out  = args.out
num  = args.num
xmax = args.xmax
ymax = args.ymax
qmax = args.qmax
fnm  = args.fnm

fnb  = '.block0.out' + out + '.'
fne = '.vtk'

for n in range(int(num) + 1):
    if args.v:
      print("Processing",n,"complete")
    vis_temp(n,prob + fnb,fne,float(qmax),int(xmax),int(ymax),fnm)

input('FINISHED. Press Enter to exit...')
