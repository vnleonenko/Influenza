#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SI-model of Baroyan-Rvachev origin with discrete time

Version history:
    * v2 - using numpy.optimize instead of iterating through k
    * v3 - including I0 into fitted params
    * v4 - considering bias in t_peak (iterating through possible value range)
    * v5 - refactor
"""

import itertools
import time
from functools import partial

import numpy as np

from core import datetime_functions as dtf
from core.methods import BaroyanSLSQPOptimizer, BaroyanLBFGSBOptimizer, BaroyanTNCOptimizer, BaroyanNelderMeadOptimizer
from core.models.baroyan_rvachev import BaroyanParams
from core.utils import get_flu_data, remove_background_incidence, get_city_name, get_filename_list, \
    get_population
from draw_data import plot_fit

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "5.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"


class Params(BaroyanParams):
    N = 3000  # epid duration
    T = 8  # disease duration for a single person
    SIZE = 25
    DISABLE_RANDOM = True
    K_RANGE = (1.02, 1.6)
    I0_RANGE = (10000.0, 10000.0)  # (0.1, 100)
    TPEAK_BIAS_RANGE = range(-7, 7)  # (-3, 3)

    # Uncomment to run BaroyanGeneticOptimizer
    # SIZE = 1
    POPULATION_SIZE = 500
    CX_PROBABILITY = 0  # 0.5
    MUT_PROBABILITY = 0.2  # 0.2
    GENERATIONS_COUNT = 50


def fit(args, population, city_mark):
    file, optimizer_cls = args
    dates, y_real = get_flu_data(file)
    data = remove_background_incidence(y_real)

    t = time.time()
    optimizer = optimizer_cls(data, population[file[-12:-8]], Params())
    optimizer.fit_one_outbreak()
    y_model, R_square_opt, k_opt, I0_opt, tpeak_bias_opt, delta = optimizer.get_results()
    elapsed_time = time.time() - t

    date_int = int(file[-12:-8] + file[-8:-6] + file[-6:-4])  # извлекается из имени файлов
    filepath = 'out/%s/' % city_mark
    filename_out_txt = 'K_out_%s_%s.txt' % (date_int, optimizer_cls.__name__)

    with open(filepath + filename_out_txt, 'ab') as f_handle:
        # f_handle.write((str(y_model) + '\n').encode())
        # f_handle.write(('Opt calc: ' + str(R_square_opt) + '\n').encode())
        np.savetxt(f_handle,
                   np.column_stack((date_int, R_square_opt, k_opt, I0_opt, tpeak_bias_opt, delta)),
                   fmt="%d %f %f %f %d %d")
        f_handle.write(('Elapsed time: %d seconds\n\n' % elapsed_time).encode())

    out_png = filepath + 'fig4_{0}_{1}_{2}.png'.format(date_int, city_mark, optimizer_cls.__name__)
    dates_new = [dtf.convertFloatToDate(x) for x in dates]
    city_name = get_city_name(city_mark)

    plot_fit(data, y_model, delta, out_png, dates_new, R_square_opt, city_name, optimizer_cls.__name__)

    return file[-12:-8]


def fit_safe(*args, **kwargs):
    try:
        return fit(*args, **kwargs)
    except Exception as e:
        return e


def invoke(files, optimizers, population, city_mark, parallel=True, safe=True):
    if safe:
        function = fit_safe
    else:
        function = fit

    if parallel:
        calc = partial(function, population=population, city_mark=city_mark)

        import multiprocessing
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        pool.map(calc, itertools.product(files, optimizers))
        pool.close()
        pool.join()

    else:
        for args in itertools.product(files, optimizers):
            function(args, population, city_mark)


def main():
    for city_mark in ['msk']:  # for three cities ,'msk','nsk']
        population = get_population(r'input_population/population_%s.csv' % city_mark)
        all_files = get_filename_list(r'FLU_rjnamm_rev/FLU_%s/' % city_mark)

        # parse_and_plot_results(city_mark, [BaroyanGeneticOptimizer], all_files)
        # invoke(all_files, [BaroyanGeneticOptimizer], population, city_mark)
        invoke(all_files,
               [BaroyanSLSQPOptimizer, BaroyanLBFGSBOptimizer, BaroyanTNCOptimizer, BaroyanNelderMeadOptimizer],
               population, city_mark)

if __name__ == '__main__':
    t0 = time.time()
    main()
    print('Total elapsed: %d seconds' % (time.time() - t0))
