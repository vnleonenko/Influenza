#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Auxiliary visualization procedures"""

import datetime

import matplotlib
import matplotlib.dates as plt_dates
import pylab as plt

import datetime_functions as dtf
from utils import get_flu_data, remove_background_incidence, get_city_name


def cut_zero_data(y_model):
    """Finds the time moment to start model data plotting"""
    i = 0
    while y_model[i] < 10 and i < len(y_model) - 1:
        i += 1
    return i


# noinspection PyPep8Naming
def plot_fit(y_real, y_model, delta, filename, flu_dates, R2, city_name, method_name):
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

    plt.title('{0}, {1} to {2} ({3})'.format(
        city_name,
        dtf.convertDateToStringMY(date_first + datetime.timedelta(days=model_beg_index)),
        dtf.convertDateToStringMY(date_first + datetime.timedelta(delta + len(y_real))),
        method_name))
    plt.grid()

    plt.savefig(filename, dpi=150, bbox_inches='tight')
    # plt.show()
    plt.close()


# noinspection PyPep8Naming
def parse_and_plot_results(city_mark, methods, all_files):
    import ast
    for method in methods:
        for file in all_files:

            dates, y_real = get_flu_data(file)
            data = remove_background_incidence(y_real)

            date_int = int(file[-12:-8] + file[-8:-6] + file[-6:-4])
            filepath = 'out25/%s/' % city_mark
            filename_out_txt = 'K_out_%s_%s.txt' % (date_int, method.__name__)
            out_txt = filepath + filename_out_txt

            with open(out_txt, 'r+') as f:
                y_model = ast.literal_eval(f.readline()[:-1])
                f.readline()
                _, R_square_opt, _, _, _, delta = f.readline()[:-1].split(' ')
                R_square_opt = float(R_square_opt)
                delta = int(delta)

                out_png = filepath + 'fig4_{0}_{1}_{2}.png'.format(date_int, city_mark, method.__name__)
                dates_new = [dtf.convertFloatToDate(x) for x in dates]
                city_name = get_city_name(city_mark)

                plot_fit(data, y_model, delta, out_png, dates_new, R_square_opt, city_name, method.__name__)
