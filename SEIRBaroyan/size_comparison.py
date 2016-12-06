#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Наша текущая задача в общей формулировке:
выяснить, можно ли добиться ускорения по сравнению с BFGS, не проиграв
в качестве калибровки (или проиграв незначительно).


    * Подвигать значения параметров алгоритмов оптимизации (options) и
        посмотреть, приведёт ли это к ускорению вычислений без потери
        качества подгонки.

    Например, здесь есть набор параметров для SLSQP (таблица внизу):
    http://www.pyopt.org/reference/optimizers.slsqp.html#optimizer-options

    * По указанной выше ссылке нашёл интересный комментарий к параметрам:
        pll_type -> STR: Parallel Implementation
        (None, ‘POA’-Parallel Objective Analysis), Default = None
        Нужно посмотреть, нет ли там встроенной параллельной реализации для
        оптимизационного алгоритма, и если это так, насколько корректно
        она работает.

    * Отдельно для SLSQP можно посмотреть, как скажется на качестве подгонки
        увеличение SIZE — например, увеличить его в два раза и погонять на
        ваших старых примерах для сравнения.
"""

import itertools
import os
import time
from functools import partial

import matplotlib
import pylab as plt

from core.optimizer import FluParams
from core.utils import get_flu_data, remove_background_incidence, get_city_name, parse_csv, get_filename_list

__author__ = "Nikita Seleznev (ne.seleznev@gmail.com)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "0.1"


INCIDENCE_ROOT = r'FLU_rjnamm_rev/FLU_%s/'
POPULATION_CSV_FILE = r'input_population/population_%s.csv'
OUT_PATH = 'benchmark/SizeComparison/'
OUT_FILE = '%04d_%04d_%s_%02d.txt'  # year, rand_seed, optimizer, size

MIN_SIZE, MAX_SIZE = 1, 50
RAND_SEEDS = [42, 420, 4200]


class Params(FluParams):
    N = 3000  # epid duration
    T = 8  # disease duration for a single person
    K_RANGE = (1.02, 1.6)
    I0_RANGE = (1.0, 1.0)  # (0.1, 100)
    TPEAK_BIAS_RANGE = range(-7, 7)  # (-3, 3)

    SIZE = None  # Iterate [1; 25]
    DISABLE_RANDOM = True
    RANDOM_SEED = None  # Iterate [42, 420, 4200]


def fit(args, filename, population, city_mark):
    params, optimizer_cls = args
    dates, y_real = get_flu_data(filename)
    data = remove_background_incidence(y_real)

    year_int = int(filename[-12:-8])
    filepath = OUT_PATH  # % city_mark
    filename_out_txt = OUT_FILE % (year_int, params.RANDOM_SEED, optimizer_cls.__name__, params.SIZE)

    if os.path.exists(filepath + filename_out_txt):
        print('PASS ' + filename_out_txt)
        return

    t0 = time.time()
    optimizer = optimizer_cls(data, population[filename[-12:-8]], params)
    optimizer.fit_one_outbreak()
    y_model, R_square_opt, k_opt, I0_opt, tpeak_bias_opt, delta = optimizer.get_results()
    elapsed_time = time.time() - t0

    with open(filepath + filename_out_txt, 'ab') as f_handle:
        f_handle.write((str(R_square_opt) + '\n').encode())
        f_handle.write(('Elapsed time: %d seconds\n\n' % elapsed_time).encode())

    print(filename_out_txt)
    return filename_out_txt


def plot_comparison(city_year_rand, opt_data, times_data, smooth):
    """
    Plot graph (R^2 x SIZE) comparison for different methods
    :param city_year_rand: str, ex. "msk,1986,0420"
    :param opt_data: dict(), key is str (Optimizer name),
        value is list of R^2, None, or dirty negative values.
    :param times_data: dict(), key is str (Optimizer name),
        value is int (number of elapsed seconds).
    :param smooth: bool, True -- draw stairway-like graph,
        False -- real data.
    """
    city_mark, year, rand_seed = city_year_rand.split(',')
    rand_seed = int(rand_seed)

    fig = plt.figure(figsize=(10, 6) if smooth else (20, 12))
    matplotlib.rcParams.update({'font.size': 14})

    ax = fig.add_subplot(111)

    x_axis = [i for i in range(MIN_SIZE, MAX_SIZE + 1)]
    colors = {'c', 'm', 'y', 'k'}
    opt_color = {
        'SLSQPOptimizer': 'b',
        'LBFGSBOptimizer': 'r',
        'TNCOptimizer': 'g',
    }

    for optimizer_name in sorted(opt_data.keys()):
        r_squares = opt_data[optimizer_name]
        if optimizer_name not in opt_color:
            opt_color[optimizer_name] = colors.pop()

        if not r_squares[MIN_SIZE]:
            r_squares[MIN_SIZE] = 0
        y_axis = [max(r_squares[MIN_SIZE], 0)]
        for i in range(MIN_SIZE + 1, MAX_SIZE + 1):
            if not r_squares[i]:
                r_squares[i] = 0
            r_squares[i] = max(r_squares[i], 0)
            if smooth:
                y_axis.append(max(r_squares[i], max(y_axis)))
            else:
                y_axis.append(r_squares[i])

        plt.plot(x_axis, y_axis, opt_color[optimizer_name] + "o-",
                 label='%s %.03f (%d sec)' % (optimizer_name, max(y_axis), times_data[optimizer_name]),
                 linewidth=2.0)

    # plt.figtext(0.15, 0.8, "$R^2 = %.3f$" % R2, fontsize=27)

    plt.xlabel('Initial parameters ' + ('number' if smooth else 'index'))
    plt.ylabel('R^2')

    plt.title('{0}, {1}. Random seed = {2}'.format(get_city_name(city_mark), year, rand_seed))
    plt.legend(loc='lower right', numpoints=1,
               prop={'size': 16}, fancybox=True, shadow=True)
    plt.grid()

    plt.savefig(OUT_PATH + 'graphs/' + ('smooth_' if smooth else '') + city_year_rand + '.png',
                dpi=450, bbox_inches='tight')
    # plt.show()
    plt.close()


def draw_all(city_mark, smooth=True):
    data = dict()  # data['msk,1986,42'] = dict();
    # dict[optimizer_class] = list(), where list[i] = R^2 or None for size i
    times = dict()

    for filename in get_filename_list(OUT_PATH):  # % city_mark):
        f_name = filename.split('/')[-1].split('.')[0]
        year, rand_seed, optimizer_name, size = f_name.split('_')
        key = ','.join([city_mark, year, rand_seed])
        with open(filename) as file:
            r_square = float(file.readline()[:-1])
            time_sec = float(file.readline().split()[-2])

        if key not in data:
            data[key] = dict()
            times[key] = dict()
        if optimizer_name not in data[key]:
            data[key][optimizer_name] = [None] * (MAX_SIZE + 1)
            times[key][optimizer_name] = 0

        data[key][optimizer_name][int(size)] = r_square
        times[key][optimizer_name] += time_sec

    for key in sorted(data.keys()):
        plot_comparison(key, data[key], times[key], smooth)


def fit_safe(*args, **kwargs):
    try:
        return fit(*args, **kwargs)
    except Exception as e:
        return e


def execute(params_list, optimizers, filename, population, city_mark, parallel=True, safe=True):
    function = fit_safe if safe else fit

    if parallel:
        calc = partial(function, filename=filename, population=population, city_mark=city_mark)

        import multiprocessing
        pool = multiprocessing.Pool(multiprocessing.cpu_count() - 2)
        pool.map(calc, itertools.product(params_list, optimizers))
        pool.close()
        pool.join()

    else:
        for args in itertools.product(params_list, optimizers):
            function(args, filename, population, city_mark)


def main():
    city_mark = 'msk'  # for three cities iterate ['msk', 'spb', 'nsk']

    population = {}  # year: population
    for item in parse_csv(POPULATION_CSV_FILE % city_mark):
        population[item[0]] = float(item[1])

    for filename in get_filename_list(INCIDENCE_ROOT % city_mark):

        params_list = []
        for rand_seed in RAND_SEEDS:
            for size in range(MIN_SIZE, MAX_SIZE + 1):
                params = Params()
                params.RANDOM_SEED = rand_seed
                params.SIZE = size
                params_list.append(params)

        # execute(params_list, [SLSQPOptimizer, LBFGSBOptimizer, TNCOptimizer],
        #         filename, population, city_mark)

    draw_all(city_mark, smooth=False)

if __name__ == '__main__':
    t0 = time.time()
    main()
    print('Total elapsed: %d seconds' % (time.time() - t0))
