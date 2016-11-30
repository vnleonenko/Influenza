import operator
import numpy as np
import fit_functions as ff
import matplotlib.pyplot as plt
import matplotlib.patches as patches

ILI_level_delta = 1.4 #1.17
relpeak_for_epid = 1819 #1819 absolute difference between the epid peak and the ILI_high_lev to be considered the epidemic

#MAGIC NUMBER - suitable for SPb only!

excess_a1 = 1.1
lack_a2 = 0.9

def extractEpidStages ( incid_list, epid_markers ):
    #Finds out ILI levels and extracts epidemic peaks

    #incid_list = data_list[...,3]
    #epid_markers contains expert opinion on epid start

    high_ILI_level_start = 75
    high_ILI_level_finish = 340
    low_ILI_level_start = 0
    low_ILI_level_finish = high_ILI_level_start

    #TODO: add curve analysis to find exactly the intervals for every given season

    ydata = incid_list

    ydata_up_init = ydata[low_ILI_level_start:low_ILI_level_finish] #a part of the curve with upper horizontal ILI graph

    a1  = ff.find_regr_coefficient(ydata_up_init)
    ydata_up_init = ydata[high_ILI_level_start:high_ILI_level_finish]
    a2 = ff.find_regr_coefficient(ydata_up_init)

    epid_peaks, epid_curves = find_epid_peaks(ydata, a2, epid_markers)#a2*ILI_level_delta)

    return a1, a2, epid_peaks, epid_curves

def suitsEpidPeakCriteria (value, level):

    if int(value) > int(level):
        return True
    else:
        return False


def suitsEpidPeakCriteria (epid_peak, ILI_level):
#criteria by rel peak

    relpeak_ill = float(epid_peak - ILI_level)#/float(ILI_level)

    if relpeak_ill>relpeak_for_epid:
        return True
    else:
        return False

def suitsEpidPeakCriteriaFormal (epid_peak_index, epid_markers_list):

    return (epid_markers_list[epid_peak_index] == 1)


def stop_condition_epi_left(value, level, eps_level, prev_el):
    #stop occurs when value is below high ILI level or the value function derivative is zero near the level
    return value < level  or ( value < eps_level*level and value-prev_el >0 )

def stop_condition_epi_right(value, level, eps_level, prev_el, edge_left):
    #stop occurs when value is below high ILI level or the value function derivative is zero near the level or right side is lower than left one
    return value < level  or value < edge_left or ( value < eps_level*level and value < eps_level*edge_left and value-prev_el > 0 )

def extractPeakElements (arr, index, ground):
    "Removing peaks, both epidemic and non-epidemic, and returning modified array along with peak incid nums and their indices"
    #!!! Assuming there's one or two peaks, - in the latter case the second is higher, - otherwise the indices will be corrupted
    #new_version: checking left border by the ILI_high level and derivative sign change, right border is set equal for the same level

    epid_curve_right = []
    epid_curve_indices_right = []
    epid_curve_left = []
    epid_curve_indices_left = []

    #scanning points in backward direction
    prev_element = arr[index] +100;
    i=index-1
    while not ( stop_condition_epi_left( arr[i], ground, ILI_level_delta, prev_element ) ) and i > 0:
        epid_curve_left.append(arr[i])
        epid_curve_indices_left.append(i)
        prev_element = arr[i]
        i=i-1

    lowest_el = arr[i+1] #epid_curve_left[-1]
    arr = list(arr)

    prev_element = arr[index] +100;
    prev_slope = 0

    #scanning points in forward direction
    while not( stop_condition_epi_right( arr[index], ground, ILI_level_delta, prev_element, lowest_el) ) and index < len(arr)-1:
        #print(arr[index])
        prev_element = arr[index]
        epid_curve_right.append(arr[index])
        epid_curve_indices_right.append(index)
        index+=1



    epid_curve_left = epid_curve_left[::-1]
    epid_curve_indices_left = epid_curve_indices_left[::-1]
    epid_curve_left.extend(epid_curve_right)
    epid_curve_indices_left.extend(epid_curve_indices_right)

    if len(epid_curve_indices_left)>0:
        arr = [ arr[i] for i in range(len(arr)) if i < epid_curve_indices_left[0] ] #removing all the data to the right of epid peak (the biggest is the most right)

    return arr, np.column_stack((epid_curve_indices_left, epid_curve_left))


def find_epid_peaks ( arr, level_delta, epid_markers):
    "returns epid peaks from the high ILI period incidence data"
    #we distinguish "just peaks" (epi criteria doesn't hold) and "epidemic peaks" (epi criteria holds)

    #level_high_bound = int(level_delta) #a2 level of high ILI
    level_ILI = int(level_delta)
    max_value = level_ILI+10 #start condition
    epid_peaks = []
    epid_curves_matr = []


    while max_value>level_ILI:

        max_value_index, max_value = max(enumerate(arr), key=operator.itemgetter(1))
        arr, epid_curve_matr = extractPeakElements( arr, max_value_index, level_ILI )

        #if suitsEpidPeakCriteria(max_value, level_ILI):   #Using our proper criteria
        if suitsEpidPeakCriteriaFormal(max_value_index, epid_markers): #Using official healthcare criteria
            epid_peaks.append(max_value)
            epid_curves_matr.append(epid_curve_matr)

    return epid_peaks, epid_curves_matr

def plotEpidemicOutbreak ( plt, col_epid ):
    #for i in range(0, len(col_epid), 1):
    i=0
    while i<len(col_epid):
        if col_epid[i] == 1.0:
            x_beg = i
            while col_epid[i] == 1.0:
                i+=1
            x_end = i-1
            plt.axvspan(x_beg, x_end, facecolor='r', alpha=0.2)
        i+=1



def plotEpidCurves ( epid_curves, days_list, x_thursday ):

    for epid_curve in epid_curves:

        plt.plot(epid_curve[...,0], epid_curve[...,1], "r", linewidth=2)
        days_list_peak = days_list[[int(i) for i in epid_curve[...,0]]]

        epid_days = epid_curve[...,0]

        x_thursday_peak = [int(i) for i in epid_days if i in x_thursday]
        x_thursday_peak_indices = [i for i in range(0, len(epid_days)) if epid_days[i] in x_thursday]

        plt.plot(x_thursday_peak, epid_curve[x_thursday_peak_indices,1], "ro", linewidth=2)
        fname_curve = 'epi_' + str(int(days_list[epid_curve[0,0]])) + '.txt'
        np.savetxt(fname_curve, np.column_stack((days_list_peak, epid_curve[...,1])), fmt="%d %d")

    plt.plot([], [], "r", label='Epidemic outbreak')
    #plt.plot([], [], "r", label='Эпидемия гриппа')
    plt.legend(loc='upper left')


#section for functions connected with transitions between ILI levels

def stop_condition_hilo(value, level, eps_level, prev_el, prev_slope):
    #Right curve is stopped by the change of the curve form for the incidence decline
    return value >=level or (value>=eps_level*level and value-prev_el <=prev_slope)

def stop_condition_lohi(value, level, eps_level, prev_el):
    #Left curve is stopped by the stop of the incidence growth
    return value >=level or (value>=eps_level*level and value-prev_el <=0)

def find_lohi_transit( ydata, a1, a2 ):
    index = 0


    prev_element = 2*ydata[index]

    #seeking the beginning
    while not(ydata[index]>a1*excess_a1 and ydata[index]>prev_element or index>len(ydata)-1):
        prev_element=ydata[index]
        index+=1

    begin_index = index

    #iterating over transition curve

    trans_curve = []
    trans_curve_index = []

    prev_element = 0.5*ydata[index]

    while not( stop_condition_lohi(ydata[index], a2, lack_a2, prev_element) or index>len(ydata)-1):
        #trans_curve.append(ydata[index])
        trans_curve_index.append(int(index))
        prev_element = ydata[index]
        index+=1

    #saving the end of a1->a2

    end_index = index


    return np.array(trans_curve_index)


def find_hilo_transit( ydata, a1, a2 ):
    index = len(ydata)-1


    prev_element = ydata[index]+100
    prev_slope = 0

    #seeking the beginning
    while not(ydata[index]>a1*excess_a1 and ydata[index]>=prev_element or index<0):
        prev_element=ydata[index]
        index-=1

    end_index = index

    #iterating over transition curve

    trans_curve = []
    trans_curve_index = []

    prev_element = ydata[index]-100

    while not( stop_condition_hilo(ydata[index], a2, lack_a2, prev_element, prev_slope) or index<0):
        #trans_curve.append(ydata[index])
        trans_curve_index.append(int(index))
        prev_element = ydata[index]
        prev_slope = ydata[index]-prev_element
        index-=1

    #saving the end of a1->a2

    begin_index = index

    #return begin_index, end_index
    return np.array(trans_curve_index[::-1])


def find_transit_curves ( ydata, a1, a2 ):
    trans_curve1_index = find_lohi_transit(ydata, a1, a2)
    trans_curve2_index = find_hilo_transit(ydata, a1, a2)
    return trans_curve1_index, trans_curve2_index


def plotLevelTransitions (trans_curve1_index, trans_curve2_index, ydata, x_thursday):
    ydata = np.array(ydata)

    if (trans_curve1_index!=[]):
        x_thursday_curve1 = [int(i) for i in trans_curve1_index if i in x_thursday]
        plt.plot(trans_curve1_index, ydata[trans_curve1_index], "g", linewidth=2, label='Level transition')
        plt.plot(x_thursday_curve1, ydata[x_thursday_curve1], "go", linewidth=4)
        #plt.plot(trans_curve1_index, ydata[trans_curve1_index], "g", linewidth=1, label='Level transition')


    if (trans_curve2_index!=[]):
        x_thursday_curve2 = [int(i) for i in trans_curve2_index if i in x_thursday]
        plt.plot(trans_curve2_index, ydata[trans_curve2_index], "g", linewidth=2)
        plt.plot(x_thursday_curve2, ydata[x_thursday_curve2], "go", linewidth=4)






