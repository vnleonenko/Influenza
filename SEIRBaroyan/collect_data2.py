#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Auxiliary script for paper, DGTS conference.
    Building surface (K, tpeak_bias, R_square)
"""

import csv
from functools import partial
import time

from scipy import array

from core.methods import BaroyanSLSQPOptimizer
from core.models.baroyan_rvachev import BaroyanParams
from core.utils import get_flu_data, remove_background_incidence, get_filename_list, \
    get_population

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"

INCIDENCE_ROOT = r'FLU_rjnamm_rev/FLU_%s/'
POPULATION_CSV_FILE = r'input_population/population_%s.csv'
DGTS_ROOT = 'benchmark/dgts/'
CSV_FILE = DGTS_ROOT + '%s_%s.csv'
PNG_FILE = DGTS_ROOT + 'graph/%s_%s.png'


class Params(BaroyanParams):
    N = 3000  # epid duration
    T = 8  # disease duration for a single person
    SIZE = 5  # 25
    DISABLE_RANDOM = True
    RANDOM_SEED = 4200
    K_RANGE = (1.02, 1.6)
    I0_RANGE = (1.0, 1.0)  # (0.1, 100)
    TPEAK_BIAS_RANGE = range(-7, 7)  # (-3, 3)


def fit(file, optimizer_cls, population, city_mark):
    dates, y_real = get_flu_data(file)
    data = remove_background_incidence(y_real)

    logger_list = []

    with open(CSV_FILE % (city_mark, file[-12:-8]), "a+") as f:
        reader = csv.reader(f)
        for idx, row in enumerate(reader):
            # print(str(idx) + ": " + str(row))
            logger_list.append((float(element) for element in row))

    optimizer = optimizer_cls(data, population[file[-12:-8]], Params(), logger_list=logger_list)
    optimizer.fit_one_outbreak()
    optimizer.get_results()

    with open(CSV_FILE % (city_mark, file[-12:-8]), 'w+', newline='') as f:
        writer = csv.writer(f)
        for row in set(logger_list):
            writer.writerow(row)

    # Drawing
    data = array(logger_list)

    xs = data[:, 0]
    ys = data[:, 1]
    zs = data[:, 2]

    try:
        from mayavi import mlab

        # Define the points in 3D space
        # including color code based on Z coordinate.
        pts = mlab.points3d(xs, ys, zs, zs)

        # Triangulate based on X, Y with Delaunay 2D algorithm.
        # Save resulting triangulation.
        mesh = mlab.pipeline.delaunay2d(pts)

        # Remove the point representation from the plot
        pts.remove()

        # Draw a surface based on the triangulation
        surf = mlab.pipeline.surface(mesh)

        # Simple plot.
        mlab.xlabel("x")
        mlab.ylabel("y")
        mlab.zlabel("z")
        mlab.show()
    except ImportError:
        import sys
        print("Unable to Import mayavi. Drawing in matplotlib", file=sys.stderr)

        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
        from matplotlib import cm
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        surf = ax.plot_trisurf(xs, ys, zs, cmap=cm.jet, linewidth=0)
        fig.colorbar(surf)

        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.zaxis.set_major_locator(MaxNLocator(5))

        fig.tight_layout()

        # plt.show()  # or:
        fig.savefig(PNG_FILE % (city_mark, file[-12:-8]))


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
    city_mark = 'msk'  # for three cities ,'msk','nsk']
    population = get_population(POPULATION_CSV_FILE % city_mark)
    all_files = get_filename_list(INCIDENCE_ROOT % city_mark)

    # Filter 2000 only
    # all_files = [x for x in all_files if "2000" in x]

    # Filter first two
    all_files = all_files[:2]

    invoke(all_files, BaroyanSLSQPOptimizer,
           population, city_mark, parallel=False, safe=False)


if __name__ == '__main__':
    t0 = time.time()
    main()
    print('Total elapsed: %d seconds' % (time.time() - t0))
