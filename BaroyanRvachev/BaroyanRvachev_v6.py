#SI-model of Baroyan-Rvachev origin with discrete time
#v2 - using numpy.optimize instead of iterating through k
#v3 - including I0 into fitted params
#v4 - considering bias in t_peak (iterating through possible value range)
#v5 - prediction on partial data
#v6 - fixed algorithm, flexible number of init seeds
import matplotlib
from scipy.optimize import minimize
import numpy as np
import pylab as plt
import datetime
import datetime_functions as dtf
import matplotlib.dates as mdates
import os, itertools, fnmatch, time
import csv
from time import time

N=3000 #epid duration
T = 8 #disease duration for a single person

def readFromCsvToList(filename):
    # return list with all data from csv file, skipping the first row (headers)
    reader = csv.reader(open(filename), delimiter=';')
    next(reader)
    res_list = list(reader)
    return res_list

def generateNumSequenceFromZero(final):
#generates the sequence of the form 0, +-1, +-2, ..., +-final
    num_list = [0]

    for i in range(1,final):
        num_list.extend([i, -i])

    return np.array(num_list)


def find_residuals( data ):
    res = 0
    mean = np.mean(data)
    for i in range(0, len(data)):
        res+=pow(data[i] - mean, 2)
    return res

def GetFluData(path):
    """
    На вход - двухколоночные данные:
     - дата в форме YYYYMMDD
     - число заболевших в течение дня
    """
    data = np.loadtxt(path)
    flu_dates = data[:,0]
    flu = data[:,1]

    return (flu_dates, flu)

def cityName(city_mark):
    if city_mark=='spb':
        city_name = 'St Petersburg'
    else:
        if city_mark=='msk':
            city_name = 'Moscow'
        else:
            if city_mark=='nsk':
                city_name = 'Novosibirsk'

    return city_name

def q( day ):
#returns infectivity by the day of infection - function q in B-R
    switcher = {
        2: 1,
        3: 0.9,
        4: 0.55,
        5: 0.3,
        6: 0.15,
        7: 0.05,
    }
    return switcher.get(day, 0) #zero infectivity by default

def m( day ):
#returns the  modifier for daily attendance of sick persons to healthcare units
    switcher = {
        0: 1.1,
        1: 1.07,
        2: 1.0,
        3: 0.95,
        4: 0.9,
        5: 0.79,
        6: 0.3,
    }
    return switcher.get(day, 0) #zero infectivity by default


def refine_data_from_raw (y):
#refines daily data using the statistics on daily attendance coefficients
    y_refined = [y[i] * m(i % 7) for i in range(0,len(y))]

    #replacing missed data by -1
    arr = np.array(y_refined)
    arr[arr<0]=-1
    return list(arr)

def remove_background_incidence ( y ):
#Considering that in the lowest incidence day the disease incidence equals background
    y_min = min(y)
    return [y[i]-y_min for i in range(0,len(y))], y_min

def doesContainEpidPeak( inc_data ):
    prev_inc = 0
    for cur_inc in inc_data:
        if cur_inc<prev_inc:
            return True #We have a peak in data (considering there's no minor peaks before the big one)
        prev_inc = cur_inc

    return False #lacking data on incidence peaks

def convert_data_to_raw (y):
#converts the model daily data back to observed data
    return [y[i] / m(i % 7) for i in range(0,len(y))]

def max_elem_index( my_list ):
#returns the index of a highest incidence
    max_value = max(my_list)
    max_index = my_list.index(max_value)
    return max_index

def sum_ill (y, t):
#summing the cumulative infectivity of the infected on the moment t
    sum = 0

    for epid_day in range(0, T):
        if t-epid_day<0:
            y_cur = 0
        else:
            y_cur = y[t-epid_day]

        #sum = sum + y[t-epid_day]*q(epid_day)
        sum = sum + y_cur*q(epid_day)
    return sum

def calculate_dist_squared(x, y, delta):
#calculating the fitting coefficient r
#x is real data, y is modeled curve
    #delta is the difference between the epidemic starts in real data and modeled curve
    sum = 0
    for i in range(delta,delta+len(x)):
        if x[i-delta]>0 and y[i]>0: #do not consider absent data which is marked by -1
            sum = sum + pow(x[i-delta] - y[i], 2)

    return sum

def calculate_peak_bias(x, y):
    x_peak = max(x)
    print("Real peak height:", x_peak)
    y_peak = max(y)
    print("Model peak height:", y_peak)
    #return abs(x_peak-y_peak)
    return y_peak/x_peak

def calculate_r(x, y, delta):
#calculating the fitting coefficient r
#x is real data, y is modeled curve
#delta is the difference between the epidemic starts in real data and modeled curve
    sum1 = 0
    sum2 = 0
    for i in range(delta,delta+len(x)):
        if x[i-delta]>0 and y[i]>0: #do not consider absent data which is marked by -1
            sum1 = sum1 + x[i-delta]*y[i]
            sum2 = sum2 + pow(y[i],2)

    return float(sum1)/float(sum2)

def calculate_s(k):
#calculating the parameter s to find the initial values of alpha and rho
    sum = 0
    for tau in range(0,T):
        sum = sum + pow(k,T-tau)*q(tau)

    return float(pow(k,T+1))/float(sum)

def MakeSimulation (alpha, lam, rho, I0):

    y = np.zeros((N+1))
    x = np.zeros((N+1))

    #initial data
    x[0] = alpha*rho
    y[0]=I0

    for t in range(0,N):
        y[t+1] = lam* x[t]*sum_ill(y,t) /rho
        x[t+1] = x[t] - y[t+1]

    return y

def CutZeroData(y_model):
#Finds the time moment to start model data plotting
    i=0
    while y_model[i]<10 and i<len(y_model)-1:
        i+=1
    return i

def PlotFit(y_real_init, y_real, y_model, delta, fname, flu_dates, R2, city_mark):
#Plotting model vs real data
    fig = plt.figure(figsize=(10, 6))
    matplotlib.rcParams.update({'font.size': 14})

    model_beg_index = CutZeroData(y_model)

    ax = fig.add_subplot(111)
    max_len = max(len(y_model), model_beg_index + len(y_real))

    date_first = flu_dates[0]-datetime.timedelta(days=delta)
    date_last = date_first + datetime.timedelta(days=max_len)

    date_range = mdates.drange(date_first,
                     date_last, datetime.timedelta(days=1))

    #plt.plot(range(model_beg_index,delta+len(y_real)), y_model[model_beg_index:delta+len(y_real)], 'g--',label='Model', linewidth = 2.0)
    #plt.plot(range(delta, delta+len(y_real)),y_real,'bo', label='Data', markersize=6)

    ###DATES on Ox######
    plt.plot_date(date_range[delta : delta + len(y_real_init)], y_real_init, "wo", label='Data', markersize=6)
    plt.plot_date(date_range[delta : delta + len(y_real)], y_real, "bo", label='Available data', markersize=6)
    plt.plot_date(date_range[model_beg_index : delta + len(y_real_init)], y_model[model_beg_index:delta+len(y_real_init)], "g--", label='Model',linewidth=2.0)

    #print([dtf.convertFloatToDate(x) for x in date_range[delta : delta + len(y_real)]])

    hfmt = mdates.DateFormatter('%d.%m')
    ax.xaxis.set_major_formatter(hfmt)

    plt.legend(loc='best',fancybox=True, shadow=True)

    plt.figtext(0.15, 0.8, "$R^2 = %.3f$" % (R2),  fontsize=27)
    ####################

    plt.ylabel('Absolute ARI incidence, cases')

    plt.title(cityName(city_mark)+", "+dtf.convertDateToStringMY(date_first+ datetime.timedelta(days=model_beg_index))+" to "+dtf.convertDateToStringMY(date_first+datetime.timedelta(delta + len(y_real))))
    plt.grid()

    plt.savefig(fname, dpi=150, bbox_inches='tight')
    #plt.show()
    plt.close()

def FindModelFit(k, rho, I0, tpeak_bias):
#Launching the simulation for a given parameter value and aligning the result to model
    s = calculate_s(k)
    alpha = 1
    lam = s
    #print(s)

    y_model = MakeSimulation(alpha, lam, rho, I0)

    #Aligning output by incidence peaks
    peak_index_real = max_elem_index(list(FluOptimizer.data))
    peak_index_model = max_elem_index(list(y_model))
    #delta = peak_index_model - peak_index_real + tpeak_bias #adding possible peak moment bias
    delta = peak_index_model - tpeak_bias #adding possible peak moment bias
    #print('Delta = ', delta)

    #######################################################################
    if delta<0:
        print("Model peak index is to the left of data peak!")
        return 10e10, [], -1
        exit()

    if peak_index_model>600:
        #print('Modeled epidemics is too long: peak index', peak_index_model)
        return 10e10, [], -1
        exit()
    #######################################################################

    #Searching for the correction coefficient
    r = calculate_r(FluOptimizer.data, y_model, delta)
    if r>1:
        r = 1
    alpha = r
    lam = s/r

    y_model = [r*y_model[i] for i in range(0, len(y_model))] #adjusting model curve for the best fit

    dist2 = calculate_dist_squared(FluOptimizer.data, y_model, delta)

    #####################
    #fig = plt.figure(figsize=(10, 6))
    #matplotlib.rcParams.update({'font.size': 14})
    #model_beg_index = CutZeroData(y_model)
    #
    #ax = fig.add_subplot(111)
    #max_len = max(len(y_model), model_beg_index + len(FluOptimizer.data))
    #
    #date_first = FluOptimizer.dates[0]-datetime.timedelta(days=delta)
    #date_last = date_first + datetime.timedelta(days=max_len)
    #
    #date_range = mdates.drange(date_first,
    #                 date_last, datetime.timedelta(days=1))

    #plt.plot_date(date_range[delta : delta + len(FluOptimizer.full_ydata)], FluOptimizer.full_ydata, "wo", label='Data', markersize=6)
    #plt.plot_date(date_range[delta : delta + len(FluOptimizer.data)], FluOptimizer.data, "bo", label='Available data', markersize=6)
    #plt.plot_date(date_range[model_beg_index : delta + len(FluOptimizer.full_ydata)], y_model[model_beg_index:delta+len(FluOptimizer.full_ydata)], "g--", label='Model',linewidth=2.0)
    #plt.show()
    ######################

    return dist2, y_model, delta

def findPeakStats( allFiles ):
    peak_day = []

    for file in allFiles:
        dates, y_real = GetFluData(file)
        peak_day.append(max_elem_index(list(y_real)))

    return peak_day


class FluOptimizer:

    full_ydata = [] #full data available without cutting
    data = [] # входной массив измерений
    rho = 0

    dates = []

    dist_zero = 0
    tpeak_bias = 0

    def __init__(self, *args):
        pass

    # функция оптимизации
    @staticmethod
    def fitFunction(k):
        #print(k[0],k[1])
        dist2, y_model, delta = FindModelFit(k[0], FluOptimizer.rho, k[1], FluOptimizer.tpeak_bias)
        #print(1- dist2/FluOptimizer.dist_zero)

        return dist2
        #return peak_bias


    def fitOneOutbreak(path, sample_size, city_mark, pop_quantity):
        dates, y_real = GetFluData(path)
        y_real, level_zero = remove_background_incidence( y_real ) #assessing baseline ARI level and subtracting from data

        real_peak = max_elem_index(y_real)
        print('Data peak: ', max_elem_index(y_real))

        if real_peak<sample_size:
            print('The peak index exceeded!')
            return

        FluOptimizer.full_ydata = y_real

        if sample_size != -1: #cutting data
            dates = dates[:sample_size]
            y_real = y_real[:sample_size]

        FluOptimizer.data = y_real

        FluOptimizer.rho = pop_quantity

        res2 = find_residuals(y_real) #Finding the squares of residuals between the real data and it math expectation

        FluOptimizer.dist_zero = res2

        k_range = (1.02,1.6)
        I0_range = (0.1, 50.0)

        peak_left  = 20
        peak_right = 80 #empirical data

        #if doesContainEpidPeak( y_real ): #narrowing the search
        #    tpeak_bias_range = (-7,7)
        #else:
            #tpeak_bias_range = (-sample_size - 3 + peak_left, -sample_size +3 + peak_right ) #-7,100 # -50
        tpeak_bias_range = ( - 3 + peak_left, 3 + peak_right ) #-7,100 # -50



        params_range = (k_range,I0_range)

        #cur_opt_fit = 0
        k_opt = 0
        R_square_opt = 0
        tpeak_bias_opt = 0
        I0_opt = 0.1
        FluOptimizer.dates = [dtf.convertFloatToDate(x) for x in dates]

        # generating uniformly distributed init values for k
        size2 = 5 #25
        init_list = []
        for i in range(0, len(params_range)):
            init_list.append(np.random.uniform(params_range[i][0], params_range[i][1], size2))
        init_list=np.array(init_list)

        for j in range(len(init_list[0,:])):
            params_init = [init_list[0,j],init_list[1,j]] #k, I0
            #print(params_init)
            ###
            if R_square_opt > 0.8 and j>1:
                print('Break!')
                break
            ###
            for tpeak_bias_cur in range(tpeak_bias_range[0], tpeak_bias_range[1]):
                t0 = time()
                print('Cur_bias:', tpeak_bias_cur)
                FluOptimizer.tpeak_bias = tpeak_bias_cur
                K = minimize(FluOptimizer.fitFunction, params_init, method='L-BFGS-B', bounds=params_range)

                fun_val = K.fun #fit value
                R_square = 1- fun_val/res2
                #peak_bias = K.fun

                print('R= ', R_square)
                #print('done, dist peak= ', K.fun)

                K1 = list(K.x) #final bunch of optimal values
                k_cur=K1[0]
                I0_cur = K1[1]

                if R_square> R_square_opt:
                #if peak_bias< peak_bias_opt:
                    k_opt = k_cur
                    I0_opt = I0_cur
                    R_square_opt = R_square
                    tpeak_bias_opt = tpeak_bias_cur
                    print('R_opt: ', R_square_opt)
                    #print(peak_bias_opt)
                t1 = time()
                print('Total time spent: %f' %(t1-t0))
                ###
                if R_square_opt > 0.9 and j>1:
                    print('Break!')
                    break
                ###
            #print(k_opt, R_square)


        if R_square_opt>0:
            dist2, y_model, delta = FindModelFit(k_opt, FluOptimizer.rho, I0_opt, tpeak_bias_opt)

            peak_height_bias = calculate_peak_bias(FluOptimizer.full_ydata+level_zero, y_model+level_zero)

            #print('Opt: ', R_square_opt)
            ###
            fun_val = dist2
            R_square_opt = 1- dist2/res2
            print('Opt calc: ', 1- dist2/res2)

            myDate = (path[-12:-8],path[-8:-6],path[-6:-4])
            fpath = 'results/BaroyanRvachev/'+city_mark+'/'
            fname_out_txt = 'K_out_'+city_mark+'.txt'

            dates_new = [dtf.convertFloatToDate(x) for x in dates]
            f_handle = open(fpath+fname_out_txt, 'ab')

            myDate = (path[-12:-8],path[-8:-6],path[-6:-4])# извлекается из имени файлов
            myDate_int = int(str(myDate[0])+str(myDate[1])+str(myDate[2]))
            print('Delta = ', delta)
            print('T_peak = ', tpeak_bias_opt)
            np.savetxt(f_handle, np.column_stack((myDate_int, sample_size, R_square_opt, tpeak_bias_opt, peak_height_bias, k_opt, I0_opt, delta, real_peak)), fmt="%d %d %f %d %f %f %f %d %d")
            f_handle.close()

            fname_out_plot = 'fig3_{0}{1}{2}_'.format(myDate[0],myDate[1],myDate[2]) + city_mark+'_'+str(sample_size)+'.png'
            PlotFit(FluOptimizer.full_ydata + level_zero, FluOptimizer.data + level_zero, y_model + level_zero, delta, fpath+fname_out_plot, dates_new, R_square_opt, city_mark)



for city_mark in ['spb']: #for three cities ,'msk','nsk']
    print(city_mark)
    population = {} # year: population
    population_list = readFromCsvToList('../../dta/population/'+city_mark+'.csv')

    for item in population_list:
            population[item[0]] = float(item[1])

    #root = r'FLU\\spb\\'
    root = '../data/incidence/'+city_mark+'/'
    allFiles =  list( itertools.chain(* [ [os.path.join(x[0],  f) for f in fnmatch.filter( x[2],"*.txt")] for x in os.walk(root) ]) )

    #peakStats = findPeakStats(allFiles)
    #print(peakStats)
    #print(min(peakStats))
    #print(max(peakStats))

    #init_data_point_numbers = [-1]
    #init_data_point_numbers = [5] #[5,10,15,20,25,30,35,40,45]#range(5,45)
    l = list(range(5,45))
    init_data_point_numbers = [58,46]#[elem for elem in l if not elem%5 ==0 ]

    for file in allFiles:
        for i in init_data_point_numbers:
            print("Init points = ", i)
            myYear = (file[-12:-8])
            #print(population[myYear])
            FluOptimizer.fitOneOutbreak(file, i, city_mark, population[myYear])


