#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Fit and draw graphs using SLSQP Baroyan-Rvachev Model implementation
"""

import datetime
import json
import os
import time
from functools import partial

import matplotlib
import matplotlib.dates as m_dates
import pylab as plt

from common import get_incidence_filenames, get_flu_data, RESULTS_PATH, get_city_name
from common import get_population
from core import datetime_functions as dtf
from core.methods import BaroyanSLSQPOptimizer
from core.models.baroyan_rvachev import BaroyanParams
from core.models.helpers import max_elem_index
from draw_data import cut_zero_data

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2017, ITMO University"
__version__ = "5.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"

CITIES = ['spb', ]  # ['msk', 'nsk', 'spb']
OUT_PATH = RESULTS_PATH + '/old/BaroyanRvachev/'
OUT_FILE = '%s_%d_%diter.json'  # city, year_month_day, size
OUT_GRAPH_FILE = 'fig5_%s_%d_%diter.pdf'  # city, year_month_day, size


class Params(BaroyanParams):
    N = 3000  # epid duration
    T = 8  # disease duration for a single person
    K_RANGE = (1.02, 1.6)
    I0_RANGE = (1.0, 1.0)  # (0.1, 100)
    TPEAK_BIAS_RANGE = range(-7, 7)  # (-3, 3)

    SIZE = 50


# noinspection PyPep8Naming
def fit(file, optimizer_cls, population, city_mark):
    dates, initial_data = get_flu_data(file)
    # initial_data = remove_background_incidence(initial_data)

    if not len(dates):  # zero output
        return

    real_peak = max_elem_index(initial_data)
    # print('Data peak: ', real_peak)

    sample_size = None  # FIXME I'm dummy, put me in args
    if sample_size and real_peak < sample_size:
        print('The peak index exceeded!')
        return

    if sample_size:  # cutting data
        dates = dates[:sample_size]
        data = initial_data[:sample_size]
    else:
        data = initial_data[:]

    date_begin = dtf.convert_float_to_date(dates[0])
    # print("date:begin", date_begin)
    date_end = dtf.convert_float_to_date(dates[-1])
    # print("date:end", date_end)

    date_int = int(file[-12:-8] + file[-8:-6] + file[-6:-4])  # извлекается из имени файлов
    filepath = OUT_PATH  # % city_mark
    filename_out_txt = OUT_FILE % (city_mark, date_int, Params.SIZE)

    # Cache results in pretty json
    results = dict()
    if os.path.exists(filepath + filename_out_txt):
        with open(filepath + filename_out_txt) as file_handle:
            results = json.load(file_handle)
            print('Using cached results ' + filename_out_txt)

    else:
        t = time.time()
        optimizer = optimizer_cls(
            data=data,
            population_quantity=population[file[-12:-8]],
            params=Params())
        optimizer.fit_one_outbreak()
        results = optimizer.get_results()
        elapsed_time = time.time() - t
        results.update({'elapsed_time': elapsed_time})

        with open(filepath + filename_out_txt, 'w') as file_handle:
            file_handle.write(
                json.dumps(
                    results, sort_keys=True,
                    indent=4, separators=(',', ': ')
                )
            )

    y_model, R_square_opt, k_opt, I0_opt,\
        tpeak_bias_opt, delta\
        = (results[key] for key in [
            "y_model", "R_square_opt", "k_opt",
            "I0_opt", "tpeak_bias_opt", "delta"])

    R_square_opt = float(R_square_opt)
    delta = int(delta)

    out_pdf = filepath + OUT_GRAPH_FILE % (city_mark, date_int, Params.SIZE)
    dates_new = [dtf.convert_float_to_date(x) for x in dates]
    city_name = get_city_name(city_mark)

    plot_fit(data, y_model, delta, out_pdf, dates_new, R_square_opt, city_name)
    print("Drawn " + out_pdf +
          "\nCalculated in " + str(int(results['elapsed_time'])) + " seconds\n")


def plot_fit(y_real, y_model, delta, filename, flu_dates, R2, city_name):
    """Plotting model vs real data"""
    fig = plt.figure(figsize=(10, 6))
    matplotlib.rcParams.update({'font.size': 14})

    model_beg_index = cut_zero_data(y_model)

    ax = fig.add_subplot(111)
    max_len = max(len(y_model), model_beg_index + len(y_real))

    date_first = flu_dates[0] - datetime.timedelta(days=delta)
    date_last = date_first + datetime.timedelta(days=max_len)

    date_range = m_dates.drange(date_first, date_last, datetime.timedelta(days=1))

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

    formatter = m_dates.DateFormatter('%d.%m')
    ax.xaxis.set_major_formatter(formatter)

    plt.figtext(0.15, 0.8, "$R^2 = %.3f$" % R2, fontsize=27)

    plt.legend(loc='best', fancybox=True, shadow=True)

    plt.ylabel('Absolute ARI incidence, cases')

    plt.title('{0}, {1} to {2}'.format(
        city_name,
        dtf.convertDateToStringMY(date_first + datetime.timedelta(days=model_beg_index)),
        dtf.convertDateToStringMY(date_first + datetime.timedelta(delta + len(y_real)))))
    plt.grid()

    plt.savefig(filename, dpi=450, bbox_inches='tight', format='pdf')
    # plt.show()
    plt.close()


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
    for city_mark in CITIES:
        population = get_population(city_mark)
        all_files = get_incidence_filenames(city_mark)

        invoke(all_files, BaroyanSLSQPOptimizer, population, city_mark, parallel=False, safe=False)

if __name__ == '__main__':
    t0 = time.time()
    main()
    print('Total elapsed: %d seconds' % (time.time() - t0))
