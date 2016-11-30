from scipy.optimize import curve_fit
import numpy as np

def func(x, a):
    return 0*x+a

def func2(x,a,b):
    return a*x+b

def sum_squares( x, y ):
    #calculating sum of squares of difference between the elements of two arrays OF EQUAL LENGTH
    sq_ar = 0;
    for i in range(0,len(x)):
        sq_ar = sq_ar + (x[i]-y[i])*(x[i]-y[i])
    return sq_ar

def delete_max_dist ( arr1, arr2 ):

    #deleting the element from arr1 which is max distant from a

    arr1 = list(arr1)

    max_dist = abs(arr1[0]-arr2[0])
    max_dist_ind = 0
    for i in range (1,len(arr1)):
        if abs(arr1[i] - arr2[i]) > max_dist:
            max_dist = abs(arr1[i] - arr2[i])
            max_dist_ind = i

    del arr1[max_dist_ind]
    return arr1

def find_regr_coefficient ( ydata ):
    "Approximate array with a horizontal line y=a2 by means of iterative least squares method with far elements removal"
    alpha = 20
    dist_sqr = 200

    sz = len(ydata)

    xdata = np.linspace(0, sz-1, sz) #range(0, sz)

    xdata_up = xdata
    ydata_up = ydata

    a2 = 0

    while dist_sqr>alpha:
        popt, pcov = curve_fit(func, xdata_up, ydata_up)
        a2  = popt
        dist_sqr = sum_squares(ydata_up, func(xdata_up,a2) )
        ydata_up = delete_max_dist(ydata_up, func(xdata_up,a2))
        xdata_up=np.linspace(0, len(ydata_up)-1, len(ydata_up))

    return a2

