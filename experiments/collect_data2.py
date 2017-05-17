#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Auxiliary script for paper, DGTS conference.
    Building surface (K, tpeak_bias, R_square)
"""

import csv
from functools import partial
import itertools
import os
from pathlib import Path
import time

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy import array

from common import get_flu_data, get_population, get_incidence_filenames, RESULTS_PATH, remove_background_incidence
from core.methods import BaroyanSLSQPOptimizer, BaroyanLBFGSBOptimizer,\
    BaroyanTNCOptimizer
from core.models.baroyan_rvachev import BaroyanParams

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"

CSV_FILE = RESULTS_PATH + '/dgts/%s/%s_%s_%s.csv'
PNG_FILE = RESULTS_PATH + '/dgts/%s/%s_%s_%s.pdf'


class Params(BaroyanParams):
    N = 3000  # epid duration
    T = 8  # disease duration for a single person
    SIZE = 25
    DISABLE_RANDOM = True
    RANDOM_SEED = 42
    K_RANGE = (1.02, 1.4)
    I0_RANGE = (1.0, 1.0)  # (0.1, 100)
    TPEAK_BIAS_RANGE = range(-7, 7)  # (-3, 3)


def fit(args, population, city_mark):
    file, optimizer_cls, params = args
    dates, y_real = get_flu_data(file)
    data = remove_background_incidence(y_real)

    logger_list = []
    params_ = (file[-12:-8], city_mark, optimizer_cls.__name__, params.SIZE)
    csv_filename = CSV_FILE % params_
    png_filename = PNG_FILE % params_

    os.makedirs(os.path.dirname(csv_filename), exist_ok=True)

    csv_file = Path(csv_filename)
    # don't calculate again the same; cutoff for rerun script
    with open(csv_filename, "r" if csv_file.is_file() else "a+") as f:
        reader = csv.reader(f)
        for idx, row in enumerate(reader):
            # print(str(idx) + ": " + str(row))
            logger_list.append(tuple(float(element) for element in row))

    if len(logger_list) == 0:  # results doesn't exists
        optimizer = optimizer_cls(data, population[file[-12:-8]], params,
                                  logger_list=logger_list)
        optimizer.fit_one_outbreak()
        optimizer.get_results()

        # store only the best R^2 values
        data_ = dict()  # (k) -> (I_0, R^2)
        for item in logger_list:
            if item[0] in data_:
                data_[item[0]].append(item[1:])
            else:
                data_[item[0]] = [item[1:]]
        logger_list = []
        for key, value in data_.items():
            I_0, R2 = max(value, key=lambda v: v[1])  # with greatest R^2
            logger_list.append((key, I_0, R2, ))

        # Save results to file
        with open(csv_filename, 'w+', newline='') as f:
            writer = csv.writer(f)
            for row in sorted(logger_list, key=lambda x: x[0]):
                writer.writerow(row)

    # Draw the graph
    data = array(logger_list)
    xs = data[:, 0]
    ys = data[:, 1]
    zs = data[:, 2]

    # try:
    #     from mayavi import mlab
    #
    #     # Define the points in 3D space
    #     # including color code based on Z coordinate.
    #     pts = mlab.points3d(xs, ys, zs, zs)
    #
    #     # Triangulate based on X, Y with Delaunay 2D algorithm.
    #     # Save resulting triangulation.
    #     mesh = mlab.pipeline.delaunay2d(pts)
    #
    #     # Remove the point representation from the plot
    #     pts.remove()
    #
    #     # Draw a surface based on the triangulation
    #     surf = mlab.pipeline.surface(mesh, extent=(0, 1, 0, 1, 0, 1))
    #
    #     # Simple plot.
    #     mlab.xlabel("K")
    #     mlab.ylabel("shift")
    #     mlab.zlabel("R^2")
    #     mlab.show()  # TODO disable
    #     mlab.savefig(png_filename[:-3] + "png")#, dpi=480)
    #
    # except ImportError:
    # import sys
    # print("Unable to Import mayavi. Drawing in matplotlib", file=sys.stderr)

    fig = plt.figure(figsize=(12, 9))
    matplotlib.rcParams.update({'font.size': 26})
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_trisurf(xs, ys, zs, cmap=cm.jet, linewidth=0)
    # fig.colorbar(surf)

    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.zaxis.set_major_locator(MaxNLocator(4))

    ax.set_xlabel('\n\nk')
    ax.set_ylabel('\n\n$\Delta_p$')
    ax.set_zlabel('\n\n$R^2$')
    plt.figtext(0.15, 0.8, "$R^2 = %.3f$" % max(zs), fontsize=32)

    fig.tight_layout()

    # plt.show()  # or:
    fig.savefig(png_filename, dpi=480, format='pdf', bbox_inches='tight',
                pad_inches=0)


def fit_safe(*args, **kwargs):
    try:
        return fit(*args, **kwargs)
    except Exception as e:
        return e


def invoke(files, optimizers, params_list, population, city_mark,
           parallel=True, safe=True):
    if safe:
        function = fit_safe
    else:
        function = fit

    if parallel:
        calc = partial(function, population=population, city_mark=city_mark)

        import multiprocessing
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        pool.map(calc, itertools.product(files, optimizers, params_list))
        pool.close()
        pool.join()

    else:
        for args in itertools.product(files, optimizers, params_list):
            function(args, population, city_mark)


def main():
    city_mark = 'spb'  # for three cities ,'msk','nsk']
    population = get_population(city_mark)
    all_files = get_incidence_filenames(city_mark)

    # Filter 2000 only
    all_files = [x for x in all_files if "2000" in x]

    # Filter first two
    # all_files = all_files[:2]
    params_list = []
    for size in [3, 9, 15]:
        params = Params()
        params.SIZE = size
        params_list.append(params)

    invoke(all_files,
           [BaroyanSLSQPOptimizer,
            BaroyanLBFGSBOptimizer,
            BaroyanTNCOptimizer],
           params_list,
           population, city_mark)  # , parallel=False, safe=False)


if __name__ == '__main__':
    t0 = time.time()
    main()
    print('Total elapsed: %d seconds' % (time.time() - t0))
