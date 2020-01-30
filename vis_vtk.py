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

def rndr(dat,var,fnmb,fnmn,min,max,n,qmax,xmx,ymx):

    x,y,z,q = dat

    fig = figure()
    img = imshow(q[var][0],extent=[-xmx,xmx,-ymx,ymx],vmin=0.0,vmax=qmax)
    fig.colorbar(img)
    title("Render of " + fnmb + ", t = " + str(n))
    xlabel("x (kpc)")
    ylabel("y (kpc)")
    savefig('rndr_' + fnmb + fnmn + '.png')
    close(fig)

    # return 0

def vtk_plt(rad,qrt,var,fnm,num):

    fig = figure()

    plot(rad[:-1],qrt,label='t = 0')

    xlim(0,20)
    ylim(0,1)
    xlabel('R (kpc)')
    ylabel(var + ' / ambient ')
    title('Ratio of ' + var + '/ambient ' + var + ' as a function of radius')

    savefig('plt_' + fnm + ('0' * (4 - omag(num,0))) + str(num) + '.png')
    close(fig)

def vis_vtk(var,num,fnb,fne,xmx,ymx,fnm,qmax):

    dat = vtk(fnb + ('0' * (4 - omag(num,0))) + str(num) + fne)
    x,y,z,q = dat

    rad = []
    qty = []

    i = 0
    j = 0
    k = 0

    while j < ymx:
        while i < xmx:
            qty.append(q[var][k,j,i])
            rad.append((x[i]**2 + y[j]**2)**0.5)

            i += 1
        i = 0
        j += 1

    qbk = max(qty)
    rndr(dat,var,fnm,('0' * (4 - omag(num,0))) + str(num),0,max(qty),num,qmax,xmx,ymx)

    qbn = stats.binned_statistic(rad,qty,'mean',bins=2046)[0]
    rbn = stats.binned_statistic(rad,qty,'mean',bins=2046)[1]

    qrt = qbn / qbk   # set ratio

    vtk_plt(rbn,qrt,var,fnm,num)
