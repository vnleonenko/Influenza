import csv
import fnmatch
import os

import itertools
import numpy as np


def parse_csv(filename):
    """Get list with all data from csv file, skipping the first row (headers)"""
    reader = csv.reader(open(filename), delimiter=';')
    next(reader)
    res_list = list(reader)
    return res_list


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


def remove_background_incidence(y):
    """Considering that in the lowest incidence day
    the disease incidence equals background+1"""
    y_min = min(y) - 1
    return [y[i] - y_min for i in range(0, len(y))]


def get_filename_list(path):
    return list(
        itertools.chain(*[
            [os.path.join(x[0], f) for f in fnmatch.filter(x[2], "*.txt")]
            for x in os.walk(path)
        ])
    )
