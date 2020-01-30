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

# Takes athdf data array collection at a single
# time step and a string of variable type to be
# working with and a string for a filename with
# path as needed NOT including file type (will
# save as with .png file extension.)

# Returns int for successful or failed run
# while saving a spatial render of whichever
# variable was passed at the given time step as
# a .png using the given name/path. Also will
# show the file in case you can see it.

def rndr(dat,var,fnmb,fnmn,min,max,n):

    x = dat['x1f']
    y = dat['x2f']

    q = dat[var]

    fig = figure()
    img = imshow(q[0],extent=[-30,30,-100,100],vmin=min,vmax=max)
    fig.colorbar(img)
    title("Render of " + fnmb + ", t = " + str(n))
    xlabel("x (kpc)")
    ylabel("y (kpc)")
    savefig('rndr_' + fnmb + fnmn + '.png')
    close(fig)

    # return 0

def athdf_plt(rad,qrt,var,fnm,num):

    fig = figure()

    plot(rad[:-1],qrt,label='t = 0')

    xlim(0,25)
    ylim(0,1)
    xlabel('R (kpc)')
    ylabel(var + ' / ambient ')
    title('Ratio of ' + var + '/ambient ' + var + ' as a function of radius')

    savefig('plt_' + fnm + ('0' * (4 - omag(num,0))) + str(num) + '.png')
    close(fig)

def vis_athdf(var,num,fnb,fne,xmx,ymx,fnm):

    dat = athdf(fnb + ('0' * (4 - omag(num,0))) + str(num) + fne)

    rad = []
    qty = []

    i = 0
    j = 0

    while j < ymx:
        while i < xmx:
            x = dat['x1f'][i]
            y = dat['x2f'][j]
            q = dat[var][0,j,i]

            qty.append(q)
            rad.append((x**2 + y**2)**0.5)

            i += 1
        i = 0
        j += 1

    qbk = max(qty)
    rndr(dat,var,fnm,('0' * (4 - omag(num,0))) + str(num),0,max(qty),num)

    qbn = stats.binned_statistic(rad,qty,'mean',bins=2046)[0]
    rbn = stats.binned_statistic(rad,qty,'mean',bins=2046)[1]

    qrt = qbn / qbk   # set ratio

    athdf_plt(rbn,qrt,var,fnm,num)
