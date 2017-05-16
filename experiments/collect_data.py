#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Auxiliary script for paper, DGTS conference
"""

import json
import time
from functools import partial

from common import get_flu_data, get_population, get_incidence_filenames, RESULTS_PATH
from core.models.seir import AbstractSEIROptimizer, SEIRParams

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "5.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"


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


# noinspection PyPep8Naming
def fit(file, optimizer_cls, population, city_mark):
    dates, initial_data = get_flu_data(file)
    # initial_data = remove_background_incidence(initial_data)

    if not len(dates):  # zero output
        return

    data = initial_data[:]

    optimizer = optimizer_cls(
        initial_data=initial_data,
        data=data,
        dates=dates,
        population_quantity=population[file[-12:-8]],
        params=Params())
    optimizer.fit_one_outbreak()
    results = optimizer.get_results()


def fit_safe(*args, **kwargs):
    try:
        return fit(*args, **kwargs)
    except Exception as e:
        return e


def invoke(files, optimizer_cls, population, city_mark, processes_count=1, safe=True):
    if safe:
        function = fit_safe
    else:
        function = fit

    if processes_count > 1:
        calc = partial(function, optimizer_cls=optimizer_cls, population=population, city_mark=city_mark)

        import multiprocessing
        pool = multiprocessing.Pool(processes_count)
        pool.map(calc, files)
        pool.close()
        pool.join()

    else:
        for file in files:
            function(file, optimizer_cls, population, city_mark)


def main():
    for city_mark in ['msk']:  # for three cities ,'msk','nsk']
        population = get_population(city_mark)
        all_files = get_incidence_filenames(city_mark)

        results = dict()
        with open(RESULTS_PATH + '/dgts/speedup.json') as f:
            results = json.load(f)

        for files_count in range(1, len(all_files) + 1):

            if str(files_count) not in results:
                results[str(files_count)] = dict()

            for processes_count in [1, 2, 4]:
                if str(processes_count) in results[str(files_count)]:
                    continue

                t0 = time.time()
                invoke(all_files[:files_count], SEIRSLSQPOptimizer, population, city_mark,
                       processes_count=processes_count, safe=False)
                results[str(files_count)][str(processes_count)] = time.time() - t0
                print('%d files, %d processes, %d seconds' % (files_count, processes_count, time.time() - t0))

                dump = json.dumps(results, sort_keys=True, indent=4, separators=(',', ': '))
                with open(RESULTS_PATH + '/dgts/speedup.json', 'w') as f_handle:
                    f_handle.write(dump)


if __name__ == '__main__':
    t0 = time.time()
    main()
    print('Total elapsed: %d seconds' % (time.time() - t0))
