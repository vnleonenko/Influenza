#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Common constants and procedures to all the experiments
"""
import csv
import fnmatch

import itertools
import os

import numpy as np

# Constant to build paths. Usage: os.path.join(BASE_DIR, 'data')
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_PATH = os.path.join(BASE_DIR, 'data')
INCIDENCE_ROOT = DATA_PATH + '/incidence/%s/'
POPULATION_CSV_FILE = DATA_PATH + '/population/%s.csv'

RESULTS_PATH = os.path.join(BASE_DIR, 'experiments/results')


def parse_csv(filename):
    """Get list with all data from csv file, skipping the first row (headers)
    """
    reader = csv.reader(open(filename), delimiter=';')
    next(reader)
    res_list = list(reader)
    return res_list


def get_population(city_code: str) -> dict():
    csv_file_path = POPULATION_CSV_FILE % city_code
    """return {year: population}"""
    population = {}
    for item in parse_csv(csv_file_path):
        population[item[0]] = float(item[1])
    return population


def get_incidence_filenames(city_mark: str) -> list():
    path = INCIDENCE_ROOT % city_mark
    return list(
        itertools.chain(*[
            [os.path.join(x[0], f) for f in fnmatch.filter(x[2], "*.txt")]
            for x in os.walk(path)
        ])
    )


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
