#operates with converted files for the seasons
#Adjusts the horizontal lines to levels of seasonal ILI
#Highlights and saves epi peaks and curves
#Shows uncorrected incidence as well
#Modified curve extraction function based on decreasing slope stop condition
#v3 New peak criteria based on squares under the graphs (highest column minus ILI high level)
#v4 Fancy graphs optimized for the workshop
#v5 Retrieving data on epid peak beginning and the corresponding temp and hum
#v6 Graph for a1_a2 levels with temperature in avg
#v7 Graph for peaks and a2 levels
#v8 = main Program refactoring to enhance readability
#main_rjnamm - transparent rectangles instead of dashed lines for the epidemic period
#_rjnamm_v2 - all artifacts connected with temp and humidity are removed

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from datetime import datetime
from scipy.optimize import curve_fit

import fit_functions as ff
import epid_peak_functions_rjnamm as epif
import datetime_functions as dtf

#service functions



def savePlot ( data_list, data_list_biased, a1, a2, epid_curves, transit_curve1, transit_curve2 ):
    #build and save one plot for a given dataset

    x_column = range(0,len(data_list), 1)

    days_list = data_list[...,0]
    incid_list = data_list[...,3]
    incid_list_biased = data_list_biased[...,3]
    epid_list = data_list[...,4]

    x_labels = days_list

    ydata = incid_list
    sz = len(ydata)

    xdata = np.linspace(0, sz-1, sz)


    #extractEpidStages

    mfig = plt.figure(figsize=(17, 6.0))
    matplotlib.rcParams.update({'font.size': 24})
    mfig.subplots_adjust(hspace=.3)

    fig1 = plt.subplot(111)
    #fig1 = mfig
    fig1.set_title('ARI incidence, Saint Petersburg, '+dtf.convertfDateToFancyString(data_list[0,0])+' to '+dtf.convertfDateToFancyString(data_list[-1,0]))

    xlabelsnew = []
    for i in x_labels:
        xlabelsnew.append(dtf.returnProperDayNames(i))

    plt.xticks(range(1,len(x_column)+1,1),xlabelsnew)

    #create an array to mark actual data (a thursday for every week)
    x_thursday = dtf.returnThursdayMarks(data_list)

    plt.plot(x_column, incid_list_biased, 'c', label='Under-reported data', linewidth=2)
    plt.plot(x_thursday, incid_list_biased[x_thursday], 'co')
    plt.plot(x_column, incid_list, 'b', linewidth=2)
    plt.plot(x_thursday, incid_list[x_thursday], 'bo')

    #plotting epidemic curves

    epif.plotEpidCurves ( epid_curves, days_list, x_thursday )
    epif.plotLevelTransitions( transit_curve1, transit_curve2, incid_list, x_thursday )
    #print('Higher ILI:', a2)

    plt.plot(xdata, ff.func(xdata,a1), "b--", label='Lower non-flu ARI level', linewidth=4)
    plt.plot(xdata, ff.func(xdata,a2), "r--", label='Higher non-flu ARI level', linewidth=4)
    plt.legend(loc='upper left', fontsize=18)

    dtf.plotSpecialDays(fig1, data_list[0,0], data_list[-1,0])
    epif.plotEpidemicOutbreak(fig1, epid_list)

    out_fname = str(year)+"-"+str(year+1)+".png"
    #plt.savefig(out_fname,dpi=300, bbox_inches='tight')

    plt.show()
    plt.close()

#years = range(1997,1998)
years = range(1987,2014)
#years = range(1988,1989)

data_list = np.loadtxt(r'out_dbases\\flu_dbase_spb.txt')
#data_list_biased = np.loadtxt('flu_dbase_biased.txt')
data_list_biased = np.loadtxt(r'out_dbases\\flu_dbase_spb_biased.txt')

a1_list = []
a2_list = []

day_end = 30
month_end = 6

#day_end  = 17
#month_end = 3

#########

for year in years:
    print(year)
    list_new = dtf.getDataArrayforGivenTimePeriod(data_list, datetime(year, 7, 1, 0, 0).date(), datetime(year+1, month_end, day_end,0,0).date())
    list_biased = dtf.getDataArrayforGivenTimePeriod(data_list_biased, datetime(year, 7, 1, 0, 0).date(), datetime(year+1, month_end, day_end,0,0).date())
    incids = list_new[...,3]
    epid_markers = list_new[...,4]

    a1, a2, epid_peaks, epid_curves = epif.extractEpidStages(incids, epid_markers)

    transit_curve1, transit_curve2 = epif.find_transit_curves ( incids, a1, a2 )

    if len(transit_curve1)>0:
        ILI_high_beg_index = int(transit_curve1[-1])
    else:
        ILI_high_beg_index=0

    if len(transit_curve2)>0:
        ILI_high_end_index = int(transit_curve2[0])
    else:
        ILI_high_end_index=0

    savePlot(list_new, list_biased, a1, a2, epid_curves, transit_curve1, transit_curve2)



