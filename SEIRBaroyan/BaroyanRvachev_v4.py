#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SI-model of Baroyan-Rvachev origin with discrete time

Version history:
    * v2 - using numpy.optimize instead of iterating through k
    * v3 - including I0 into fitted params
    * v4 - considering bias in t_peak (iterating through possible value range)
    * v5 - refactor
"""

import csv
import datetime
import fnmatch
import itertools
import os
import sys
from functools import partial

import matplotlib
import matplotlib.dates as plt_dates
import numpy as np
import pylab as plt
from scipy.optimize import minimize

import datetime_functions as dtf

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "5.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"

N = 3000  # epid duration
T = 8  # disease duration for a single person
K_RANGE = (1.02, 1.6)
I0_RANGE = (10000.0, 10000.0)  # (0.1, 100)
TPEAK_BIAS_RANGE = range(-7, 7)  # (-3, 3)


def parse_csv(filename):
    """Get list with all data from csv file, skipping the first row (headers)"""
    reader = csv.reader(open(filename), delimiter=';')
    next(reader)
    res_list = list(reader)
    return res_list


def find_residuals(data_list):
    """Finding the squares of residuals between the real data and it math expectation"""
    res = 0
    mean = np.mean(data_list)
    for i in range(0, len(data_list)):
        res += pow(data_list[i] - mean, 2)
    return res


def get_flu_data(path):
    """На вход - двухколоночные данные:
        - дата в форме YYYYMMDD
        - число заболевших в течение дня
    """
    data = np.loadtxt(path)
    flu_dates = data[:, 0]
    flu = data[:, 1]
    return flu_dates, flu


def get_city_name(city_mark):
    if city_mark == 'spb':
        return 'St Petersburg'
    elif city_mark == 'msk':
        return 'Moscow'
    elif city_mark == 'nsk':
        return 'Novosibirsk'
    return 'Unknown'


def q(day):
    """returns infectivity by the day of infection - function q in B-R"""
    switcher = {
        2: 1,
        3: 0.9,
        4: 0.55,
        5: 0.3,
        6: 0.15,
        7: 0.05,
    }
    return switcher.get(day, 0)  # zero infectivity by default


# def m(day):
#    """returns the  modifier for daily attendance of sick persons to healthcare units"""
#    switcher = {
#        0: 1.1,
#        1: 1.07,
#        2: 1.0,
#        3: 0.95,
#        4: 0.9,
#        5: 0.79,
#        6: 0.3,
#    }
#    return switcher.get(day, 0)  # zero infectivity by default


# def refine_data_from_raw(y):
#     """refines daily data using the statistics on daily attendance coefficients"""
#     y_refined = [y[i] * m(i % 7) for i in range(0,len(y))]
#
#     #replacing missed data by -1
#     arr = np.array(y_refined)
#     arr[arr<0]=-1
#     return list(arr)


def remove_background_incidence(y):
    """Considering that in the lowest incidence day
    the disease incidence equals background+1"""
    y_min = min(y) - 1
    return [y[i] - y_min for i in range(0, len(y))]

# def convert_data_to_raw (y):
#     """converts the model daily data back to observed data"""
#     return [y[i] / m(i % 7) for i in range(0,len(y))]


def max_elem_index(my_list):
    """returns the index of a highest incidence"""
    return my_list.index(max(my_list))


def sum_ill(y, t):
    """summing the cumulative infectivity of the infected on the moment t"""
    sum_ = 0

    for epid_day in range(0, T):
        if t - epid_day < 0:
            y_cur = 0
        else:
            y_cur = y[t - epid_day]

        # sum_ = sum_ + y[t-epid_day]*q(epid_day)
        sum_ += y_cur * q(epid_day)
    return sum_


def calculate_dist_squared(x, y, delta):
    """calculating the fitting coefficient r
    x is real data, y is modeled curve
    delta is the difference between the epidemic starts in real data and modeled curve
    """
    sum_ = 0
    for i in range(delta, delta + len(x)):
        if x[i-delta] > 0 and y[i] > 0:  # do not consider absent data which is marked by -1
            sum_ += pow(x[i-delta] - y[i], 2)
    return sum_


def calculate_peak_bias(x, y):
    x_peak = max(x)
    y_peak = max(y)
    return abs(x_peak - y_peak)


def calculate_r(x, y, delta):
    """calculating the fitting coefficient r
    x is real data, y is modeled curve
    delta is the difference between the epidemic starts in real data and modeled curve
    """
    sum1 = 0
    sum2 = 0
    for i in range(delta, delta + len(x)):
        if x[i-delta] > 0 and y[i] > 0:  # do not consider absent data which is marked by -1
            sum1 += x[i-delta] * y[i]
            sum2 += pow(y[i], 2)
    return float(sum1) / float(sum2)


def calculate_s(k):
    """calculating the parameter s to find the initial values of alpha and rho"""
    sum_ = 0
    for tau in range(0, T):
        sum_ += pow(k, T - tau) * q(tau)

    return float(pow(k, T+1)) / float(sum_)


# noinspection PyPep8Naming
def make_simulation(alpha, lam, rho, I0):
    y = np.zeros((N+1))
    x = np.zeros((N+1))

    # initial data
    x[0] = alpha * rho
    y[0] = I0

    for t in range(0, N):
        y[t+1] = lam * x[t] * sum_ill(y, t) / rho
        # print(y[t+1])
        x[t+1] = x[t] - y[t+1]
    return y


def cut_zero_data(y_model):
    """Finds the time moment to start model data plotting"""
    i = 0
    while y_model[i] < 10 and i < len(y_model)-1:
        i += 1
    return i


# noinspection PyPep8Naming
def plot_fit(y_real, y_model, delta, filename, flu_dates, R2, city_name):
    """Plotting model vs real data"""
    fig = plt.figure(figsize=(10, 6))
    matplotlib.rcParams.update({'font.size': 14})

    model_beg_index = cut_zero_data(y_model)

    ax = fig.add_subplot(111)
    max_len = max(len(y_model), model_beg_index + len(y_real))

    date_first = flu_dates[0] - datetime.timedelta(days=delta)
    date_last = date_first + datetime.timedelta(days=max_len)

    date_range = plt_dates.drange(date_first, date_last, datetime.timedelta(days=1))

    # plt.plot(range(model_beg_index, delta + len(y_real)),
    #          y_model[model_beg_index: delta+len(y_real)],
    #          'g--',label='Model', linewidth = 2.0)
    # plt.plot(range(delta, delta + len(y_real)),
    #          y_real,
    #          'bo', label='Data', markersize=6)

    # DATES on Ox
    plt.plot_date(date_range[model_beg_index: delta + len(y_real)],
                  y_model[model_beg_index: delta + len(y_real)],
                  "g--", label='Model', linewidth=2.0)
    plt.plot_date(date_range[delta: delta + len(y_real)],
                  y_real,
                  "bo", label='Data', markersize=6)

    # print([dtf.convertFloatToDate(x) for x in date_range[delta : delta + len(y_real)]])

    formatter = plt_dates.DateFormatter('%d.%m')
    ax.xaxis.set_major_formatter(formatter)

    plt.legend(loc='best', fancybox=True, shadow=True)

    plt.figtext(0.15, 0.8, "$R^2 = %.3f$" % R2, fontsize=27)

    plt.ylabel('Absolute ARI incidence, cases')

    plt.title('{0}, {1} to {2}'.format(
        city_name,
        dtf.convertDateToStringMY(date_first + datetime.timedelta(days=model_beg_index)),
        dtf.convertDateToStringMY(date_first + datetime.timedelta(delta + len(y_real)))))
    plt.grid()

    plt.savefig(filename, dpi=150, bbox_inches='tight')
    # plt.show()
    plt.close()


# noinspection PyPep8Naming
def find_model_fit(k, rho, I0, tpeak_bias, data):
    """Launching the simulation for a given parameter value and aligning the result to model"""
    s = calculate_s(k)
    alpha = 1
    lam = s
    # print(s)

    y_model = make_simulation(alpha, lam, rho, I0)

    # Aligning output by incidence peaks
    peak_index_real = max_elem_index(list(data))
    peak_index_model = max_elem_index(list(y_model))
    delta = peak_index_model - peak_index_real + tpeak_bias  # adding possible peak moment bias

    #######################################################################
    if delta < 0:
        sys.stderr.write("Model peak index is to the left of data peak!\n")
        return 10e10, [], -1, -1
    #######################################################################

    # Searching for the correction coefficient
    r = calculate_r(data, y_model, delta)
    # alpha = r
    # lam = s/r

    y_model = [r * y_item for y_item in y_model]  # adjusting model curve for the best fit
    # print(y_model)

    dist2 = calculate_dist_squared(data, y_model, delta)

    peak_bias = calculate_peak_bias(data, y_model)

    # len_data = len(FluOptimizer.data)
    # plt.plot(range(0, len_data), FluOptimizer.data)
    # plt.plot(range(0, len_data), y_model[:len_data])
    # plt.show()

    return dist2, y_model, delta, peak_bias


# noinspection PyPep8Naming
def fit_function(params, rho, tpeak_bias, data):
    k, I0 = params
    dist2, y_model, delta, peak_bias = find_model_fit(k, rho, I0, tpeak_bias, data)
    return dist2


# noinspection PyPep8Naming
class FluOptimizer:
    def __init__(self, data, population_quantity):
        self.data = data
        self.rho = population_quantity
        self.res2 = find_residuals(self.data)

        self.k_opt = None
        self.I0_opt = None
        self.R_square_opt = 0
        self.tpeak_bias_opt = None

    def fit_one_outbreak(self):
        params_range = (K_RANGE, I0_RANGE)

        # k_opt, R_square_opt, I0_opt, tpeak_bias_opt = 0, 0, 0, 0

        # generating unifromly distributed init values for k
        size2 = 2
        np.random.seed(42)
        init_list = [np.random.uniform(param[0], param[1], size2) for param in params_range]
        init_list = np.array(init_list)

        for j in range(len(init_list[0, :])):
            params_init = [init_list[0, j], init_list[1, j]]  # k, I0
            print(params_init)

            for tpeak_bias_cur in TPEAK_BIAS_RANGE:
                print(tpeak_bias_cur)
                partial_function = partial(fit_function, rho=self.rho, tpeak_bias=tpeak_bias_cur, data=self.data)
                K = minimize(partial_function, params_init, method='L-BFGS-B', bounds=params_range)

                fun_val = K.fun  # fit value
                R_square = 1 - fun_val/self.res2
                print('done, R= ', R_square)

                # if peak_bias < peak_bias_opt:
                if R_square > self.R_square_opt:
                    k_cur, I0_cur = list(K.x)  # final bunch of optimal values
                    self.k_opt = k_cur
                    self.I0_opt = I0_cur
                    self.R_square_opt = R_square
                    self.tpeak_bias_opt = tpeak_bias_cur
                    # print(R_square_opt)
                    # print(peak_bias_opt)
            # print(k_opt, R_square)

    def save_results(self, dates, date_int, out_txt, out_png, city_name):
        dist2, y_model, delta, peak_bias = find_model_fit(
            self.k_opt, self.rho, self.I0_opt, self.tpeak_bias_opt, self.data)

        print(y_model)
        # print('Opt: ', R_square_opt)
        # fun_val = dist2
        R_square_opt = 1 - dist2/self.res2
        print('Opt calc: ', R_square_opt)

        with open(out_txt, 'ab') as f_handle:
            np.savetxt(f_handle,
                       np.column_stack((date_int, R_square_opt, self.k_opt, self.I0_opt, self.tpeak_bias_opt, delta)),
                       fmt="%d %f %f %f %d %d")

        dates_new = [dtf.convertFloatToDate(x) for x in dates]
        plot_fit(self.data, y_model, delta, out_png, dates_new, R_square_opt, city_name)


def main():
    for city_mark in ['msk']:  # for three cities ,'msk','nsk']
        population = {}  # year: population
        population_list = parse_csv(r'input_population/population_%s.csv' % city_mark)

        for item in population_list:
            population[item[0]] = float(item[1])

        root = r'FLU_rjnamm_rev/FLU_%s/' % city_mark
        all_files = list(
            itertools.chain(*[
                [os.path.join(x[0], f) for f in fnmatch.filter(x[2], "*.txt")]
                for x in os.walk(root)
            ])
        )

        # init_data_point_numbers = [-1]
        # init_data_point_numbers = range(5, 45)

        for file in all_files:
            # for i in init_data_point_numbers:
            #     print("Init points = ", i)
            my_year = (file[-12:-8])
            # print(population[my_year])
            dates, y_real = get_flu_data(file)
            data = remove_background_incidence(y_real)

            optimizer = FluOptimizer(data, population[my_year])
            optimizer.fit_one_outbreak()

            filepath = 'out25/%s/' % city_mark
            filename_out_txt = 'K_out_%s.txt' % city_mark
            out_txt = filepath + filename_out_txt
            date_int = int(file[-12:-8] + file[-8:-6] + file[-6:-4])  # извлекается из имени файлов
            out_png = 'fig3_{0}_{1}.png'.format(date_int, city_mark)
            city_name = get_city_name(city_mark)
            optimizer.save_results(dates, date_int, out_txt, out_png, city_name)

if __name__ == '__main__':
    main()
