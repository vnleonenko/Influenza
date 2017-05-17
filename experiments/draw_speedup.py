#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Auxiliary script for paper, DGTS conference.
    Draw speedup graph
"""

import json
import time

import matplotlib
import matplotlib.pyplot as plt

from common import RESULTS_PATH

SPEEDUP_FILE = RESULTS_PATH + '/dgts/speedup.json'
OUTPUT_FILE = RESULTS_PATH + '/dgts/speedup.pdf'


def main():
    data = {}
    with open(SPEEDUP_FILE) as f:
        data = json.load(f)

    speedups = dict()
    for size, measurements in data.items():
        if int(size) == 1:
            continue  # pass trivial case

        one_process = float(measurements["1"])
        for process_count, seconds in measurements.items():
            if int(process_count) == 1:
                continue  # speedup for 1 process === 1.0

            try:
                speedups[int(process_count)][int(size)] = one_process / float(seconds)
            except KeyError:
                speedups[int(process_count)] = {int(size): one_process / float(seconds)}

    fig = plt.figure(figsize=(10, 6))  # if smooth else (20, 12))
    matplotlib.rcParams.update({'font.size': 20})

    ax = fig.add_subplot(111)

    sizes = next(iter(speedups.values())).keys()
    x_axis = [i for i in range(min(sizes), max(sizes) + 1)]
    colors = {'c', 'm', 'y', 'k'}
    opt_color = {
        2: 'b', 4: 'r', 8: 'g',
    }

    for process_count, measurements in speedups.items():
        speedup_list = [measurements[key] for key in sorted(measurements.keys())]
        if process_count not in opt_color:
            opt_color[process_count] = colors.pop()

        plt.plot(x_axis, speedup_list, opt_color[process_count] + "o-",
                 label='%d processes speedup' % (process_count),
                 linewidth=2.0)

    plt.xlabel('Time periods')
    plt.ylabel('Speedup')

    plt.legend(loc='lower right', numpoints=1,
               prop={'size': 16}, fancybox=True, shadow=True)
    plt.grid()

    plt.savefig(OUTPUT_FILE, dpi=450, format='pdf', bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    t0 = time.time()
    main()
    print('Total elapsed: %d seconds' % (time.time() - t0))
