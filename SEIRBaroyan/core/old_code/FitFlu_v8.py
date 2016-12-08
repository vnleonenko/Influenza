# -*- coding: cp1251 -*-
"""
 - Моделирование процесса распространения гриппа на основе системы ДУ
 - Идентификация параметров модели на основе данных измерений градиентным методом локального поиска
 - v2 - corrected version (starting from the first measurement point if it's to the left of the model curve start
 - v3 - batch version to check the predictive force of the model
 - v4 - changed algorithm for fitting model curve to data - peak adjustment replaced by the discrete moving of the epid modeling starting moment
 - v5 - modified version to distinguish incomplete data without information about the peaks
 - v5b - starting level (coef = 0.9) is moved to fit function
 - v5cb - working on improved predictive ability
 - v6 - finding out the height and the date of an outbreak + fancy pics, rjnamm
 - v6b - still rjnamm, drawing S, E, I, R populations on the graph without baseline addition
 - v6_rjnamm - corrections according to the RJNAMM article review (new infected calc mechanism via E, R2 calculation routine)
 - v6_rjnamm_aN - experimental version using non-immune ratio instead of vert scaling and peak level to fit the curve
 - v8 - initial parameter variation added (like in _aN_intervals)
"""
# ________ Flu Model ________
# S - Susceptible (восприимчивые к болезни)
# E - Exposed (латентные - первые заразившиеся)
# I - Infected (инфицированные с признаками болезни)
# R - Recovered (выздоровевшие)
# ---- ФОРМА УРАВНЕНИЙ ----
# dS/dt  = -k1*S*I              <-- текущее здоровых
# dE/dt =  +k1*S*I - k2*E       <-- текущее инфицированных латентных
# dI/dt  = +k2*E - k3*I         <-- текущее инфицированны
# dR/dt  = +k3*I                <--- ВЫЗДОРОВЕЛО ВСЕГО (ТОЛЬКО РОСТ):
# ---------- СМЕРТНОСТЬ ----------
# dR/dt  = +k4*T                <--- УМЕРЛО ВСЕГО (ТОЛЬКО РОСТ) -- не используется

import pylab as p
import numpy as np
import datetime
import csv
import matplotlib
from matplotlib import rc
from itertools import cycle
from DrawData import *
import datetime_functions as dtf
import matplotlib.dates as mdates
import sys,os, itertools,shutil,fnmatch,time
from time import time

lines = ["-","--","-.",":"]
linecycler = cycle(lines)
font = {'family': 'Verdana','weight': 'normal'}
rc('font', **font)

k1,k2,k3 = [0.6, 0.81, 0.5 ]

min_cases_show = 10 #minimum incidence value to show to the right in the model trajectory (value for crop)

MIN_EPID_LEN = 100
MAX_EPID_LEN = 150
EPID_LEN_PENALTY = 2000000

# k1 = 0.5    # S -> E
# k2 = 0.25   # E -> I
# k3 = 0.3    # I -> T

def readFromCsvToList(filename):
    # return list with all data from csv file
    reader = csv.reader(open(filename), delimiter=';')
    next(reader)
    res_list = list(reader)
    return res_list

def readPopulation (city_mark):
    population = {} # year: population
    population_list = readFromCsvToList('population_msk.csv')

    for item in population_list:
        population[item[0]] = float(item[1])

def modelEpidemicDuration( M_data ):
#finds the duration of epidemic from the data curve

    i=0
    while M_data[i]*20000<1 and i<len(M_data)-1:
        i+=1

    epid_beg = i

    i = len(M_data)-1

    while M_data[i]*20000<1 and i>0:
        i-=1

    epid_end = i

    duration = max(0, epid_end - epid_beg)
    #print(duration)

    return duration

# --- данные (измерения) ---
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

# ---- данные ----
def ShowFluStat(figure, fluData_init, flu_dates, DAY,M, data_left, data_right, size, K_opt, optimal_value, R2, city_mark):
    #p.figure(figsize=(19.2, 12.0))
    #p.figure()

    ###Calculating the first and the last date of the range
    ax = figure.add_subplot(111)
    max_len = max(len(M), data_left + len(fluData_init))
    print(flu_dates)
    print(len(flu_dates))

    date_first = flu_dates[0]
    date_last = date_first + datetime.timedelta(days=max_len)

    ##Generating data range
    hfmt = mdates.DateFormatter('%d.%m')
    ax.xaxis.set_major_formatter(hfmt)
    date_range = mdates.drange(date_first,
                     date_last, datetime.timedelta(days=1))

    print("Model data ", len(M))

    if(M is not None):
        p.plot_date(date_range[:len(M)], M, "g--", label='Model',linewidth=2.0)
        #p.plot_date(date_range[:len(M)], S, "c", label='Susceptible',linewidth=4.0)
        #p.plot_date(date_range[:len(M)], E, "r", label='Exposed',linewidth=4.0)
        #p.plot_date(date_range[:len(M)], I, "m", label='Infected',linewidth=4.0)
        #p.plot_date(date_range[:len(M)], R, "y", label='Recovered',linewidth=4.0)

    print("Init data: ", len(fluData_init))
    print("Range: ", len(date_range[data_left : data_left + len(fluData_init)]))


    p.plot_date(date_range[data_left : data_left + len(fluData_init)], fluData_init, "o", color = 'lightgray', label='Full dataset',markersize=5)
    p.plot_date(date_range[data_left: data_right], DAY[data_left:data_right], "bo", label='Data', markersize=5)

    print('Data peak: ', max(DAY))
    print('Data peak date: ', mdates.num2date(date_range[list(DAY).index(max(DAY))]))


    print('Model peak: ', max(M))
    print('Model peak date: ', mdates.num2date(date_range[list(M).index(max(M))]))
    #print(mdates.num2date(date_range))

    p.legend(loc='best',fancybox=True, shadow=True)
    #p.text(date_range[5], max(fluData_init), "k1= {0}, k2= {1}, k3 = {2}".format(K_opt[1],K_opt[2],K_opt[3]))
    #p.text(date_range[5], max(fluData_init), "k1= %.2f, k2= %.2f, k3 = %.2f, fit = %.2E" % (K_opt[1],K_opt[2],K_opt[3], int(optimal_value)))
    p.figtext(0.15, 0.8, "$R^2 = %.3f$" % (R2),  fontsize=27)

    flu_dates = [dtf.returnSomeDaysNameFromDate(x) for x in flu_dates]

    #p.ylabel('Absolute ARI incidence, cases')

    #p.xlabel('days')

    if city_mark=='spb':
        city_name = 'Saint Petersburg'
    else:
        if city_mark=='msk':
            city_name = 'Moscow'
        else:
            if city_mark=='nsk':
                city_name = 'Novosibirsk'


    p.title(city_name + ", "+dtf.convertDateToStringMY(date_first)+" to "+dtf.convertDateToStringMY(date_last))
    p.ylabel('ARI incidence')
    p.grid()

def max_elem_index( my_list ):
#returns the index of a highest incidence
    my_list = list(my_list)
    max_value = max(my_list)
    max_index = my_list.index(max_value)
    return max_index

def calculate_peak_bias(x, y):
    return max(y)/max(x)

def calculate_tpeak_bias(x, y):
    peak_index_real = max_elem_index(x)
    peak_index_model = max_elem_index(y)
    delta = peak_index_model - peak_index_real
    return delta

#def doesContainEpidPeak( inc_data ):
#    prev_inc = 0
#    for cur_inc in inc_data:
#        if cur_inc<prev_inc:
#            return True #We have a peak in data (considering there's no minor peaks before the big one)
#        prev_inc = cur_inc
#
#    return False #lacking data on incidence peaks

# ________ общая процедура моделирования ________
def MakeSimulation(K,N): 
    # ________ ДУ ________
    inInf,k1,k2,k3 = K
    #N = N*3 # удвоим интервал моделирования

    t =  np.linspace(0, N-1, N)
    def dX_dt(X, t=0): # Return parameters of the infected populations.
        S,E,I,T = X 
        return np.array([-k1*S*I,
                      +k1*S*I - k2*E,
                      +k2*E - k3*I,
                      +k3*I])
    # ________ интегрируем ДУ ________
    from scipy import integrate    
    #inInf = 0.01    # начальное число зараженных
    # переделать на zeros
    X0 = np.array([1.-inInf, inInf, 0.0,0.0])  # initials conditions: S,E,I,T
    X,infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
    S,E,I,T = X.T

    return toInf(S,E), S[1:], E[1:], I[1:], T[1:]

# S - результат моделирования
def toInf(S,E):
    INF = [-S[i+1]+S[i] - E[i+1]+E[i] for i in range(len(S)-1)] # число зараженых за период (new cases), модель
    return INF

# --- отрезать лишние значения и произвести перенормировку---
def shiftModelToData(DAY, M, shift, s_ratio, pop_total):
#Shifting the model curve and cutting the points which are outside the data boundaries

    #Calculating the prospected model epid length under this parameter set
    param = modelEpidemicDuration(M)
    #print("Epi_duration: ", param)

    shiftM1 = shift             # столько отрезаем от начала модели (сдвиг модели влево)
    shiftM2 = len(DAY)+shift    # столько отрезаем от конца модели

    #print("Shifts: ", shiftM1, shiftM2, shiftM2 - shiftM1)

    M_new = M[shiftM1:shiftM2]
    #print(len(M))
    #print(len(M_new))

    #M_new = M_new * s_ratio * float(pop_total)
    M_new = [i* s_ratio * float(pop_total) for i in M_new]

    return (DAY,M_new, param)

# the same procedure modified for the sake of visualisation (no model trajectories cut, hence the output arrays have different sizes)---
def shiftDataToModel(dates,DAY,M, shift, s_ratio,  S, E, I, R, pop_total):
#Moving data points to the right for the sake of plotting model vs data
#Also scaling model incidence along with E and I
    N_zeros = len(M) - len(DAY) #considering the case when model data is always of greater length than the real data

    date_begin = dtf.convertFloatToDate(dates[0])
    date_end = dtf.convertFloatToDate(dates[-1])
    dates_new = []
    data_left = shift

    DAY_new = []

    DAY_new.extend([0]*shift)
    DAY_new.extend(DAY)
    data_right = len(DAY_new)
    DAY_new = list(DAY_new)

    DAY_new.extend([0]*(N_zeros-shift))

    date_begin_adj = date_begin - datetime.timedelta(days=int(shift))
    date_end_adj = date_end + datetime.timedelta(days=int(N_zeros-shift))
    M_new = M

    daysbetween1 = date_begin - date_begin_adj
    dates_range1 = [date_begin_adj + datetime.timedelta(days=x) for x in range(0, daysbetween1.days)]

    daysbetween2 = date_end_adj - date_end
    dates_range2 = [date_end + datetime.timedelta(days=x) for x in range(0, daysbetween2.days)]

    dates_new.extend(dates_range1)
    dates_new.extend([dtf.convertFloatToDate(x) for x in dates])
    dates_new.extend(dates_range2)

    # Calculating the old parameter scaling_coef for the comparison purposes
    maxMeas = max(DAY)                          # максимальное значение (измерения)
    maxModel = max(M)                           # максимальное значение (модель)
    scaling_coef = s_ratio*pop_total*maxModel / maxMeas
    #print(scaling_coef)

    #M_new = M_new * s_ratio * pop_total
    M_new = [i* s_ratio * float(pop_total) for i in M_new]

    #S_scaled = S* s_ratio * pop_total
    #E_scaled = E* s_ratio * pop_total
    #I_scaled = I* s_ratio * pop_total
    #R_scaled = R* s_ratio * pop_total

    #removing the redundant zeros to the right of M_new

    current_end = len(M_new)-1

    while M_new[current_end] < min_cases_show and current_end>data_right: #the value of M_new is less than a fixed number
        current_end-=1


    return (dates_new[:current_end+1], DAY_new[:current_end+1],M_new[:current_end+1]), data_left, data_right, scaling_coef


def find_residuals( data ):
    res = 0
    mean = np.mean(data)
    for i in range(0, len(data)):
        res+=pow(data[i] - mean, 2)
    return res


from scipy.optimize import minimize


class FluOptimizer:

    init_data = [] #full array of incidence points to measure the model predictive force
    data = [] # входной массив измерений - (v5b) хранится в "не обрезанном" по уровню виде, обрезается в fit-функции
    level_zero = 0
    population_total = 0 #quantity of population for a given season
    N_model = 0 #number of simulation points
    start_shift = 0
    cur_optimal_value = -1

    date_begin_min = (datetime.datetime.now()).date()
    date_end_max = datetime.date(1800, 1, 1)

    def __init__(self, *args):
        pass
            
    # функция оптимизации
    @staticmethod
    def fitFunction(K):
        #print(K)
        data = FluOptimizer.data
        N =  FluOptimizer.N_model #len(data)
        pop_total = FluOptimizer.population_total
        INF, S, E, I, R = MakeSimulation(K[:4],N)                   # процедура моделирования (incidence)

        s_ratio = K[4]
        shifting_coef = K[5]

        shift = FluOptimizer.level_zero*shifting_coef

        data = data - shift

        #print(scaling_coef)
        ModData = shiftModelToData(data,INF, FluOptimizer.start_shift, s_ratio, pop_total)
        DAY,M, modelEpidemicDuration = ModData[0],ModData[1],ModData[2]

        fitDif1 = sum(pow(DAY-M,2))

        return fitDif1
    
    def fitOneOutbreak(path, sample_size, shift_range, city_mark, pop_quantity):

        fd = GetFluData(path) #Taking all the data from the file

        if fd[0] == []: #zero output
            return

        dates, fluData = fd[0], fd[1]
        FluOptimizer.population_total = pop_quantity

        real_peak = max_elem_index(fluData)
        print('Data peak: ', real_peak)

        if real_peak<sample_size:
            print('The peak index exceeded!')
            return

        #FluOptimizer.level_zero = min(fluData)

        FluOptimizer.N_model = len(fluData)*3
        FluOptimizer.init_data = fluData

        if sample_size != -1: #cutting data
            dates = dates[:sample_size]
            fluData = fluData[:sample_size]

        FluOptimizer.level_zero = min(fluData)

        #print(dates[0])
        date_begin = dtf.convertFloatToDate(dates[0])
        print("date:begin",date_begin)
        date_end = dtf.convertFloatToDate(dates[-1])
        print("date:end",date_end)
        #ShowFluStat(fluData)

        FluOptimizer.data = fluData

        #shifting_coef_range = (0.8, 1.0) #0.999
        s_ratio = (0.001, 1.0) #ratio of susceptible from the total population
        scaling_coef_opt = 0
        cur_optimal_shift = 0
        #if sample_size == -1 or doesContainEpidPeak(fluData): #Analysing full epi data or the data contains a peak

        #bnds = ((0.0001, 0.00011),(0.5, 0.8), (0.4, 0.9),(0.3, 0.9), scaling_range, shifting_coef_range) #rjnamm_rev
        #bnds = ((0.0001, 0.00011),(0.000001, 10.0), (0.000001, 10.0),(0.000001, 10.0), scaling_range, shifting_coef_range) #rjnamm_rev2_ext
        #bnds = ((0.0001, 0.0001001),(0.0000001, 50.0), (0.0000001, 50.0),(0.0000001, 50.0), scaling_range, shifting_coef_range) #rjnamm_rev2_ext_big
        #bnds = ((0.0001, 0.00011),(0.095, 0.75), (0.2, 8.0),(0.08, 0.34), scaling_range, shifting_coef_range) #rjnamm_rev2 real_biol

        #bnds = ((0.0001, 0.0001001),(0.0000001, 50.0), (0.0000001, 50.0),(0.0000001, 50.0), s_ratio, shifting_coef_range) #rjnamm_rev2_ext_big
        #bnds = ((0.0001, 0.00011),(0.0000001, 50.0), (0.2, 8.0),(0.08, 0.34), s_ratio, shifting_coef_range) #rjnamm_rev2 real_biol

        #raro
        shifting_coef_range = (0.7, 1.0) #0.999
        bnds = ((0.0001, 0.0001001),(0.0000001, 50.0), (0.39, 0.390001),(0.133, 0.1330001), s_ratio, shifting_coef_range) #rjnamm_rev2 real_biol


        #initK = [0.00045801916356709254, 0.606337408461603, 0.49311295274678785, 0.25945012150235969, 1.0, 0.95]
        #initK = [0.0001, 0.606337408461603, 0.49311295274678785, 0.25945012150235969, 1.0, 1.0]
        #initK = [0.0001, 0.606337408461603, 0.49311295274678785, 0.36, 0.5, 1.0] #5a, 6a

        res2 = find_residuals(fluData) #Finding the squares of residuals between the real data and it math expectation

        #generating unifromly distributed init values for s_ratio and k1
        size2 = 5 # 25
        R_square_opt = -1
        init_list = []
        init_list.append(np.random.uniform(0.0000001, 50.0, size2)) #k1
        #init_list.append(np.random.uniform(0.2, 8.0, size2)) #k2
        #init_list.append(np.random.uniform(0.08, 0.34, size2)) #k3“
        init_list.append(np.random.uniform(0.001, 1.0, size2))  #s_ratio
        init_list.append(np.random.uniform(0.7, 1.0, size2)) #shift_coef

        init_list = np.array(init_list)

        ModData = []
        FluOptimizer.cur_optimal_value = -1

        #print(init_list[0,:])

        for j in range(len(init_list[0,:])):
            #initK = [0.0001, init_list[0,j], init_list[1,j], init_list[2,j], init_list[3,j], 0.8]
            initK = [0.0001, init_list[0,j], 0.39, 0.133, init_list[1,j], init_list[2,j]]
            print("J:",j)
            ###
            if R_square_opt>0.8 and j>1:
                print('Break!')
                break

            ###

            for cur_shift in shift_range:
                t0 = time()
                #print(cur_shift)
                FluOptimizer.start_shift = cur_shift
                K = minimize(FluOptimizer.fitFunction,initK,method='L-BFGS-B', bounds=bnds  ) # SLSQP | L-BFGS-B ,bounds=bnds
                K1 = list(K.x) #final bunch of optimal values
                fun_val = K.fun #fit value
                R_square = 1- fun_val/res2

                print("Cur R: {0} Opt R: {1}".format(R_square, R_square_opt))

                if (fun_val< FluOptimizer.cur_optimal_value or FluOptimizer.cur_optimal_value == -1) and R_square>0:
                    #print("Set")
                    #print(fun_val)
                    s_ratio_opt = K1[4]
                    K_opt = K1[:4]
                    shifting_coef_cur = K1[5]
                    cur_optimal_shift = cur_shift
                    FluOptimizer.cur_optimal_value = fun_val
                    R_square_opt = R_square
                    INF, S, E, I, R = MakeSimulation(K_opt,FluOptimizer.N_model)               # launching simulation with optimal param set for visualization purpose
                    ModData, data_left, data_right, scaling_coef_opt = shiftDataToModel(dates, fluData - FluOptimizer.level_zero*shifting_coef_cur, INF, cur_shift, s_ratio_opt, S, E, I, R, FluOptimizer.population_total)

                    ###
                    if R_square_opt > 0.9 and j>1:
                        print('Break!')
                        break
                    ###

                t1 = time()
                print('Total time spent: %f' %(t1-t0))

        # print("Value: ", FluOptimizer.cur_optimal_value)
        # print("Shift: ", cur_optimal_shift)
        # print("K_opt: ", K_opt)
        # print("Vertical shift coef: ", shifting_coef_cur)
        # print("Scaling coef: ", scaling_coef_opt)

        if R_square_opt>0: #saving only the results that make sense
            vert_shift = FluOptimizer.level_zero*shifting_coef_cur
            peak_bias = calculate_peak_bias(FluOptimizer.init_data+vert_shift,ModData[2]+vert_shift)
            tpeak_bias = calculate_tpeak_bias(FluOptimizer.init_data+vert_shift,ModData[2]+vert_shift)

            fpath = 'out\\'+city_mark+'\\'
            fname_out_txt = 'K_out_'+city_mark+'.txt'
            f_handle = open(fpath+fname_out_txt, 'ab')
            myDate = (path[-12:-8],path[-8:-6],path[-6:-4])# извлекается из имени файлов

            myDate_str = str(myDate[0])+"/"+str(myDate[1])+"/"+str(myDate[2])
            myDate_int = int(str(myDate[0])+str(myDate[1])+str(myDate[2]))
            np.savetxt(f_handle, np.column_stack((myDate_int, sample_size, R_square_opt, tpeak_bias, peak_bias, cur_optimal_shift, K_opt[0], K_opt[1], K_opt[2], K_opt[3], shifting_coef_cur, s_ratio_opt, scaling_coef_opt )), fmt="%d %d %f %d %f %d %f %f %f %f %f %f %f")
            f_handle.close()

            #p.figure()
            figure = p.figure(figsize=(10, 6))
            matplotlib.rcParams.update({'font.size': 18})

            ShowFluStat(figure, FluOptimizer.init_data, ModData[0], ModData[1]+vert_shift,ModData[2]+vert_shift, data_left, data_right, sample_size, K_opt, FluOptimizer.cur_optimal_value, R_square_opt, city_mark)
            #fname = "fig3_{0}{1}{2}_nsk.pdf".format(myDate[0],myDate[1],myDate[2])
            fname_out_plot = 'fig3_{0}{1}{2}_'.format(myDate[0],myDate[1],myDate[2]) + city_mark+'_'+str(sample_size)+'.png'
            p.savefig(fpath+fname_out_plot, dpi=150, bbox_inches='tight')
            #p.show()  # --- вывод данных на экран
            p.close()




for city_mark in ['nsk']:
    population = {} # year: population
    population_list = readFromCsvToList(r'input_population\\population_'+city_mark+'.csv')

    for item in population_list:
            population[item[0]] = float(item[1])

    #root = r'FLU\\spb\\'
    root = r'FLU_rjnamm_rev\\FLU_'+city_mark+'\\'
    allFiles =  list( itertools.chain(* [ [os.path.join(x[0],  f) for f in fnmatch.filter( x[2],"*.txt")] for x in os.walk(root) ]) )


    #try:
    #    os.remove('K_out_aN.txt')
    #except OSError:
    #    pass


    init_data_point_numbers = [-1]
    #init_data_point_numbers = range(60,68)

    #shift_range = range(2,60) extended
    shift_range = range(5,54)
    #shift_range = range(10,14)

    #for file in allFiles:
    #    for i in init_data_point_numbers:

    for i in init_data_point_numbers:
        for file in allFiles:
            print("Init points = ", i)
            myYear = (file[-12:-8])
            #print(population[myYear])
            FluOptimizer.fitOneOutbreak(file, i, shift_range, city_mark, population[myYear])




