#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SEIR Model implementation for benchmark purposes
"""

import datetime
import itertools
import time
from functools import partial

import matplotlib
import matplotlib.dates as m_dates
import numpy as np
import pylab

from core import datetime_functions as dtf
from core.models.helpers import max_elem_index
from core.models.seir import AbstractSEIROptimizer, SEIRParams
from core.utils import get_flu_data, get_city_name, get_filename_list, \
    get_population

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "5.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"

INCIDENCE_ROOT = r'FLU_rjnamm_rev/FLU_%s/'
POPULATION_CSV_FILE = r'input_population/population_%s.csv'
OUT_PATH = 'benchmark/SEIRModel/'
OUT_FILE = '%d_%02d.txt'  # year_month_day, size


class Params(SEIRParams):
    SIZE = 5
    DISABLE_RANDOM = True
    RANDOM_SEED = 42

    SHIFT_RANGE = (5, 54)  # unused
    K0_RANGE = (0.0001, 0.0001001)  # ???
    K1_RANGE = (0.0000001, 50.0)  # ???
    K2_RANGE = (0.39, 0.390001)  # ???
    K3_RANGE = (0.133, 0.1330001)  # ???
    S_RATIO_RANGE = (0.001, 1.0)  # ratio of susceptible from the total population
    SHIFT_COEFFICIENT_RANGE = (0.7, 1.0)


class SEIRSLSQPOptimizer(AbstractSEIROptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        from scipy.optimize import minimize
        result = minimize(function, minimize_params, method='L-BFGS-B', bounds=minimize_params_range)
        return result.fun, tuple(result.x)


def plot_flu_statistics(initial_data, flu_dates, DAY, M, data_left, data_right,
                        sample_size, r_square, city_mark, date_int):
    figure = pylab.figure(figsize=(10, 6))
    matplotlib.rcParams.update({'font.size': 18})
    # p.figure(figsize=(19.2, 12.0))
    # p.figure()

    # Calculating the first and the last date of the range
    ax = figure.add_subplot(111)
    max_len = max(len(M), data_left + len(initial_data))
    print(flu_dates)
    print(len(flu_dates))

    date_first = flu_dates[0]
    date_last = date_first + datetime.timedelta(days=max_len)

    # Generating data range
    formatter = m_dates.DateFormatter('%d.%m')
    ax.xaxis.set_major_formatter(formatter)
    date_range = m_dates.drange(date_first, date_last, datetime.timedelta(days=1))

    print("Model data ", len(M))

    if M is not None:
        pylab.plot_date(date_range[:len(M)], M, "g--", label='Model', linewidth=2.0)
        #p.plot_date(date_range[:len(M)], S, "c", label='Susceptible',linewidth=4.0)
        #p.plot_date(date_range[:len(M)], E, "r", label='Exposed',linewidth=4.0)
        #p.plot_date(date_range[:len(M)], I, "m", label='Infected',linewidth=4.0)
        #p.plot_date(date_range[:len(M)], R, "y", label='Recovered',linewidth=4.0)

    print("Init data: ", len(initial_data))
    print("Range: ", len(date_range[data_left : data_left + len(initial_data)]))

    pylab.plot_date(date_range[data_left: data_left + len(initial_data)],
                    initial_data, "o", color='lightgray',
                    label='Full dataset', markersize=5)
    pylab.plot_date(date_range[data_left: data_right],
                    DAY[data_left:data_right], "bo",
                    label='Data', markersize=5)

    print('Data peak: ', max(DAY))
    print('Data peak date: ', m_dates.num2date(date_range[list(DAY).index(max(DAY))]))

    print('Model peak: ', max(M))
    print('Model peak date: ', m_dates.num2date(date_range[list(M).index(max(M))]))
    # print(mdates.num2date(date_range))

    pylab.legend(loc='best', fancybox=True, shadow=True)
    # p.text(date_range[5], max(fluData_init), "k1= {0}, k2= {1}, k3 = {2}".format(K_opt[1],K_opt[2],K_opt[3]))
    # p.text(date_range[5], max(fluData_init),
    #   "k1= %.2f, k2= %.2f, k3 = %.2f, fit = %.2E" % (K_opt[1],K_opt[2],K_opt[3], int(optimal_value)))
    pylab.figtext(0.15, 0.8, "$R^2 = %.3f$" % r_square, fontsize=27)

    flu_dates = [dtf.returnSomeDaysNameFromDate(x) for x in flu_dates]

    # p.ylabel('Absolute ARI incidence, cases')
    # p.xlabel('days')

    city_name = get_city_name(city_mark)

    pylab.title('{}, {} to {}'.format(city_name, dtf.convertDateToStringMY(date_first),
                                      dtf.convertDateToStringMY(date_last)))
    pylab.ylabel('ARI incidence')
    pylab.grid()

    #####################################################33

    # fname = "fig3_{0}{1}{2}_nsk.pdf".format(myDate[0],myDate[1],myDate[2])
    filename_out_plot = '{0}_{1}_{2}.png'.format(city_mark, date_int, (sample_size or ''))
    pylab.savefig(OUT_PATH + 'graphs/' + filename_out_plot, dpi=150, bbox_inches='tight')
    # p.show()  # --- вывод данных на экран
    pylab.close()


# noinspection PyPep8Naming
def fit(file, optimizer_cls, population, city_mark):
    dates, initial_data = get_flu_data(file)
    # initial_data = remove_background_incidence(initial_data)

    if not len(dates):  # zero output
        return

    real_peak = max_elem_index(initial_data)
    print('Data peak: ', real_peak)

    sample_size = None  # FIXME I'm dummy, put me in args
    if sample_size and real_peak < sample_size:
        print('The peak index exceeded!')
        return

    if sample_size:  # cutting data
        dates = dates[:sample_size]
        data = initial_data[:sample_size]
    else:
        data = initial_data[:]

    date_begin = dtf.convertFloatToDate(dates[0])
    print("date:begin", date_begin)
    date_end = dtf.convertFloatToDate(dates[-1])
    print("date:end", date_end)

    t = time.time()
    optimizer = optimizer_cls(
        initial_data=initial_data,
        data=data,
        dates=dates,
        population_quantity=population[file[-12:-8]],
        params=Params())
    optimizer.fit_one_outbreak()
    results = optimizer.get_results()
    y_model, R_square_opt, k_opt, s_ratio_opt, shift_coefficient_opt,\
        shift_opt, S, E, I, R, peak_bias, tpeak_bias,\
        model_data, data_left, data_right, scaling_coefficient\
        = (results[key] for key in [
            "y_model", "R_square_opt", "k_opt",
            "s_ratio_opt", "shift_coefficient_opt",
            "shift_opt", "S", "E", "I", "R", "peak_bias", "tpeak_bias",
            "model_data", "data_left", "data_right", "scaling_coefficient"])
    elapsed_time = time.time() - t

    date_int = int(file[-12:-8] + file[-8:-6] + file[-6:-4])  # извлекается из имени файлов
    filepath = OUT_PATH  # % city_mark
    filename_out_txt = OUT_FILE % (date_int, Params().SIZE)
    with open(filepath + filename_out_txt, 'ab') as f_handle:
        np.savetxt(f_handle, np.column_stack((date_int, sample_size or -1, R_square_opt, tpeak_bias, peak_bias,
                                              shift_opt, k_opt[0], k_opt[1], k_opt[2], k_opt[3],
                                              shift_coefficient_opt, s_ratio_opt, scaling_coefficient)),
                   fmt="%d %d %f %d %f %d %f %f %f %f %f %f %f")
        f_handle.write(('Elapsed time: %d seconds\n\n' % elapsed_time).encode())

    level_zero = min(data)
    vertical_shift = level_zero * shift_coefficient_opt
    plot_flu_statistics(initial_data, model_data[0],
                        model_data[1] + vertical_shift, model_data[2] + vertical_shift,
                        data_left, data_right,
                        sample_size, R_square_opt,
                        city_mark, date_int)


def fit_safe(*args, **kwargs):
    try:
        return fit(*args, **kwargs)
    except Exception as e:
        return e


def invoke(files, optimizer_cls, population, city_mark, parallel=True, safe=True):
    if safe:
        function = fit_safe
    else:
        function = fit

    if parallel:
        calc = partial(function, optimizer_cls=optimizer_cls, population=population, city_mark=city_mark)

        import multiprocessing
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        pool.map(calc, files)
        pool.close()
        pool.join()

    else:
        for file in files:
            function(file, optimizer_cls, population, city_mark)


def main():
    for city_mark in ['msk']:  # for three cities ,'msk','nsk']
        population = get_population(r'input_population/population_%s.csv' % city_mark)
        all_files = get_filename_list(r'FLU_rjnamm_rev/FLU_%s/' % city_mark)[:1]

        invoke(all_files, SEIRSLSQPOptimizer, population, city_mark, parallel=False, safe=False)

if __name__ == '__main__':
    t0 = time.time()
    main()
    print('Total elapsed: %d seconds' % (time.time() - t0))
