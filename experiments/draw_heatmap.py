#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Auxiliary script for paper, DGTS conference.
    Building contour map (K, tpeak_bias, R_square)
"""
import csv
import os
from functools import partial
from pathlib import Path

import pylab as plt
from scipy.interpolate import griddata
import numpy as np
import matplotlib
from astropy.convolution import convolve
from astropy.convolution.kernels import Gaussian2DKernel
import matplotlib.colors as mcolors

from common import RESULTS_PATH, get_flu_data, remove_background_incidence, get_population, get_incidence_filenames
from models.baroyan_rvachev import AbstractBaroyanOptimizer, BaroyanParams

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2017, ITMO University"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"

# Assert it has structure dgts/YYYY/city_method_iterations.csv
DGTS_RESULTS_ROOT = os.path.join(RESULTS_PATH, 'dgts')
CSV_FILE = os.path.join(DGTS_RESULTS_ROOT, 'graph', '%s_%s_%s.csv')


class Params(BaroyanParams):
    N = 3000  # epid duration
    T = 8  # disease duration for a single person
    SIZE = 250
    K_RANGE = (1.02, 1.4)
    I0_RANGE = (1.0, 1.0)
    TPEAK_BIAS_RANGE = range(-7, 7)


class DummyBaroyanOptimizer(AbstractBaroyanOptimizer):

    def optimize(self, function, minimize_params, minimize_params_range):
        return function(minimize_params), minimize_params

    def fit_one_outbreak(self):
        logger_list = []
        params_range = (self.params.K_RANGE, self.params.I0_RANGE)

        K_min, K_max = self.params.K_RANGE
        I0_min, I0_max = self.params.I0_RANGE
        size = self.params.SIZE
        init_params = zip(
            np.linspace(K_min, K_max, size),
            np.linspace(I0_min, I0_max, size)
        )

        for params in init_params:
            for tpeak_bias_cur in self.params.TPEAK_BIAS_RANGE:
                function = self._get_fit_function(tpeak_bias_cur)

                value, args = self.optimize(function, params, params_range)
                k_cur, I0_cur = args

                R_square = 1.0 - value/self.res2
                logger_list.append((k_cur, tpeak_bias_cur, max(R_square, 0.0)))
        return logger_list


def fit(file, population, city_mark):
    """
    Computes R^2 value for product of given Ks and deltas in parameters
    Creates .csv files: [(K, delta, R^2), ..., (K, delta, R^2)]
    """
    dates, initial_data = get_flu_data(file)
    data = remove_background_incidence(initial_data)

    logger_list = []
    csv_filename = CSV_FILE % (file[-12:-8], city_mark, Params().SIZE)

    os.makedirs(os.path.dirname(csv_filename), exist_ok=True)

    csv_file = Path(csv_filename)
    # don't calculate again the same; cutoff for rerun script
    with open(csv_filename, "r" if csv_file.is_file() else "a+") as f:
        reader = csv.reader(f)
        for idx, row in enumerate(reader):
            # print(str(idx) + ": " + str(row))
            logger_list.append(tuple(float(element) for element in row))

    if len(logger_list) == 0:  # results doesn't exists
        optimizer = DummyBaroyanOptimizer(
            data=data,
            population_quantity=population[file[-12:-8]],
            params=Params())
        logger_list = optimizer.fit_one_outbreak()

        # Save results to file
        with open(csv_filename, 'w+', newline='') as f:
            writer = csv.writer(f)
            for row in sorted(logger_list, key=lambda x: x[0]):
                writer.writerow(row)
    return logger_list


def draw(K, delta, R2, output_file):
    # min_K, max_K = (1.02, 1.4)
    # min_delta, max_delta = -7, 6
    # extent = [min_K, max_K, min_delta, max_delta]

    # for bins in range(15, 16):  # (20,50):
    #     for stddev in [1.1]:#, 5.0, 10.0]:  # , 4.0, 5.0]:
    fig, ax = plt.subplots(figsize=(10, 4))
    matplotlib.rcParams.update({'font.size': 20})

    # heatmap, xedges, yedges = np.histogram2d(x=K, y=delta, bins=bins, weights=R2)
    # # c = mcolors.ColorConverter().to_rgb
    # # rvb = make_colormap([c('white'), c('yellow'), 0.33, c('yellow'), c('red'), 0.90, c('red')])
    # # heatmap = convolve(heatmap.transpose()[::-1], Gaussian2DKernel(stddev=stddev))
    # plt.imshow(heatmap, extent=extent, interpolation='bilinear', alpha=0.5, aspect='auto',
    #            cmap=plt.get_cmap('jet'))

    # Convert from pandas dataframes to numpy arrays
    X, Y, Z, = np.array([]), np.array([]), np.array([])
    for k, d, r2 in zip(K, delta, R2):
        X = np.append(X, k)
        Y = np.append(Y, d)
        Z = np.append(Z, r2)

    # create x-y points to be used in heatmap
    xi = np.linspace(X.min(), X.max(), len(X))
    yi = np.linspace(Y.min(), Y.max(), len(Y))

    # Z is a matrix of x-y values
    zi = griddata((X, Y), Z, (xi[None, :], yi[:, None]), method='cubic',
                  fill_value=0.0)

    # Create the contour plot
    CS = plt.contourf(xi, yi, zi, 500, cmap=plt.cm.rainbow,
                      vmax=1.0, vmin=0.0)
    plt.colorbar()

    # # Plot points
    # plt.scatter(X, Y)

    plt.title('Contour graph for $K$, $\delta$, $R^2$')
    ax.autoscale(False)
    # cb = plt.colorbar()
    # cb.set_label('don\'t know, %')
    plt.xlabel(r'$K$', fontsize=30)
    plt.ylabel(r'$\delta$', fontsize=30)
    # plt.savefig(output_file + '_%02d_%.0f' % (bins, stddev) + '.png', dpi=150, bbox_inches='tight')
    plt.savefig(output_file + '.png', dpi=150, bbox_inches='tight')
    # plt.savefig(output_file + '.pdf', dpi=150, bbox_inches='tight')
    # plt.show()
    plt.close()


def fit_safe(*args, **kwargs):
    try:
        return fit(*args, **kwargs)
    except Exception as e:
        return e


def invoke(files, params_list, population, city_mark,
           parallel=True, safe=True):
    if safe:
        function = fit_safe
    else:
        function = fit

    if parallel:
        calc = partial(function, population=population, city_mark=city_mark)

        import multiprocessing
        pool = multiprocessing.Pool(multiprocessing.cpu_count() - 2)
        pool.map(calc, files)  # itertools.product(files, optimizers, params_list))
        pool.close()
        pool.join()

    else:
        for args in files:  # itertools.product(files, optimizers, params_list):
            function(args, population, city_mark)


def main():
    city_mark = 'spb'
    population = get_population(city_mark)
    all_files = get_incidence_filenames(city_mark)

    # Uncomment to compute 2000 only
    # all_files = [x for x in all_files if '2000' in x]

    # Uncomment to change SIZE
    # Params.SIZE = 1000

    # Uncomment to change K range
    # Params.K_RANGE = (1.02, 1.17)

    # Comment it not to launch computation again
    invoke(all_files,
           None,
           population, city_mark)  # , parallel=False, safe=False)

    # Draw all the data we can found
    for root, dirs, files in os.walk(os.path.join(DGTS_RESULTS_ROOT, 'graph')):
        for file in files:
            if file.endswith(".csv"):
                # Uncomment to re-draw only 2000 with size 1000
                # if file[:4] != '2000' or file[-8:-4] != '1000':
                #     continue

                filename = os.path.join(root, file)
                data = np.loadtxt(filename, delimiter=',')
                K, delta, R2 = data[:, 0], data[:, 1], data[:, 2]
                output_file = filename[:-4]  # cut '.csv'
                draw(K, delta, R2, output_file)

if __name__ == '__main__':
    main()
