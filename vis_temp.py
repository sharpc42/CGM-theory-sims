#!/usr/bin/env python

from athena_read import *
from numpy import *
from matplotlib.pyplot import *
from scipy import stats
from argparse import *

def omag(amt,ord):
    if amt / 10 < 1.0:
        return ord
    else:
        return omag(amt / 10,ord + 1)

# Takes vtk data array collection at a single
# time step and a string of variable type to be
# working with and a string for a filename with
# path as needed NOT including file type (will
# save as with .png file extension.)

# Returns int for successful or failed run
# while saving a spatial render of whichever
# variable was passed at the given time step as
# a .png using the given name/path. Also will
# show the file in case you can see it.

def rndr(dat,fnmb,fnmn,min,max,n,xmx,ymx):

    x,y,z,q = dat
    temp = 1.5 * q['press'] / q['rho']

    fig  = figure()
    img  = imshow(temp[0],extent=[-xmx,xmx,-ymx,ymx],vmin=min,vmax=max)

    fig.colorbar(img)
    title("Render of Temp, t = " + str(n))
    xlabel("x (kpc)")
    ylabel("y (kpc)")
    savefig('rndr_temp' + fnmn + '.png')
    close(fig)

    # return 0

def temp_plt(rad,qrt,fnm,num):

    fig = figure()

    plot(rad[:-1],qrt,label='t = 0')

    xlim(0,25)
    ylim(0,1)
    xlabel('R (kpc)')
    ylabel('temp / ambient ')
    title('Ratio of temp / ambient temp as a function of radius')

    savefig('plt_' + fnm + ('0' * (4 - omag(num,0))) + str(num) + '.png')
    close(fig)

def vis_temp(num,fnb,fne,qmx,xmx,ymx,fnm):

    dat = vtk(fnb + ('0' * (4 - omag(num,0))) + str(num) + fne)
    x,y,z,q = dat
    temp = 1.5 * q['press'] / q['rho']

    rad = []
    qty = []

    i = 0
    j = 0
    k = 0

    while j < ymx:
        while i < xmx:
            qty.append(1.5 * q['press'][k,j,i] / q['rho'][k,j,i])
            rad.append((x[i]**2 + y[j]**2)**0.5)

            i += 1
        i = 0
        j += 1

    qbk = max(qty)
    rndr(dat,fnm,('0' * (4 - omag(num,0))) + str(num),0.00001,qmx,num,xmx,ymx)

    qbn = stats.binned_statistic(rad,qty,'mean',bins=2046)[0]
    rbn = stats.binned_statistic(rad,qty,'mean',bins=2046)[1]

    qrt = qbn / qbk   # set ratio

    temp_plt(rbn,qrt,fnm,num)
