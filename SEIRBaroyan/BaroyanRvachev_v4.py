#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SI-model of Baroyan-Rvachev origin with discrete time

Version history:
    * v2 - using numpy.optimize instead of iterating through k
    * v3 - including I0 into fitted params
    * v4 - considering bias in t_peak (iterating through possible value range)
    * v5 - refactor
"""

import fnmatch
import itertools
import os
from functools import partial

import numpy as np
from scipy.optimize import minimize

import datetime_functions as dtf
from draw_data import plot_fit
from optimizer import FluParams, FluOptimizer
from utils import get_flu_data, remove_background_incidence, get_city_name, parse_csv

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "5.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"


class Params(FluParams):
    N = 3000  # epid duration
    T = 8  # disease duration for a single person
    SIZE = 25
    # DISABLE_RANDOM = True
    K_RANGE = (1.02, 1.6)
    I0_RANGE = (10000.0, 10000.0)  # (0.1, 100)
    TPEAK_BIAS_RANGE = range(-7, 7)  # (-3, 3)


class BFGSOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='BFGS', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class LBFGSBOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='L-BFGS-B', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class SLSQPOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='SLSQP', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class NelderMeadOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='Nelder-Mead', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class PowellOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='Powell', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class CGOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='CG', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class TNCOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='TNC', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class COBYLAOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='COBYLA', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


def calc_stuff(file, population, city_mark, optimizer_clazz):
    # for i in init_data_point_numbers:
    #     print("Init points = ", i)
    # print(population[my_year])
    dates, y_real = get_flu_data(file)
    data = remove_background_incidence(y_real)

    optimizer = optimizer_clazz(data, population[file[-12:-8]], Params())
    optimizer.fit_one_outbreak()
    y_model, R_square_opt, k_opt, I0_opt, tpeak_bias_opt, delta = optimizer.get_results()

    filepath = 'out25/%s/' % city_mark
    filename_out_txt = 'K_out_%s_%s.txt' % (city_mark, optimizer_clazz.__name__)
    out_txt = filepath + filename_out_txt
    date_int = int(file[-12:-8] + file[-8:-6] + file[-6:-4])  # извлекается из имени файлов

    with open(out_txt, 'ab') as f_handle:
        f_handle.write((str(y_model) + '\n').encode())
        f_handle.write(('Opt calc: ' + str(R_square_opt) + '\n').encode())
        np.savetxt(f_handle,
                   np.column_stack((date_int, R_square_opt, k_opt, I0_opt, tpeak_bias_opt, delta)),
                   fmt="%d %f %f %f %d %d")

    out_png = filepath + 'fig4_{0}_{1}_{2}.png'.format(date_int, city_mark, optimizer_clazz.__name__)
    dates_new = [dtf.convertFloatToDate(x) for x in dates]
    city_name = get_city_name(city_mark)

    plot_fit(data, y_model, delta, out_png, dates_new, R_square_opt, city_name, optimizer_clazz.__name__)

    return file[-12:-8]


def calc_stuff_safe(optimizer_clazz, file, population, city_mark):
    try:
        return calc_stuff(file, population, city_mark, optimizer_clazz)
    except Exception as e:
        return e


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

        # parse_and_plot_results(city_mark, [TNCOptimizer], all_files)
        # for method in [NelderMeadOptimizer, SLSQPOptimizer]:
        #     for file in all_files:
        #         calc_stuff(file, population, city_mark, method)

        calc_parallel = partial(calc_stuff_safe, population=population, city_mark=city_mark, file=all_files[0])

        import multiprocessing
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        pool.map(calc_parallel, [LBFGSBOptimizer, SLSQPOptimizer, TNCOptimizer])  # NelderMeadOptimizer
        pool.close()
        pool.join()

if __name__ == '__main__':
    import time
    t0 = time.time()
    main()
    print(time.time() - t0, 'seconds')
