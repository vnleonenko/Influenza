#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Operates with converted files for the seasons
Adjusts the horizontal lines to levels of seasonal ILI
Highlights and saves epi peaks and curves
Shows uncorrected incidence as well
Modified curve extraction function based on decreasing slope stop condition

Version history:
    * v3 New peak criteria based on squares under the graphs
        (highest column minus ILI high level)
    * v4 Fancy graphs optimized for the workshop
    * v5 Retrieving data on epidemic peak beginning and the corresponding
        temp and hum
    * v6 Graph for a1_a2 levels with temperature in avg
    * v7 Graph for peaks and a2 levels
    * v8 = main Program refactoring to enhance readability
    * main_rjnamm - transparent rectangles instead of dashed lines for the
        epidemic period
    * rjnamm_v2 - all artifacts connected with temp and humidity are removed
    * v11 refactor, add plot_separate() to show all cities for all years
        add plot_at_same_graph() to show multiple cities epidemic curve at
        the same graph for all years
"""

from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import matplotlib.lines as plt_lines
import matplotlib.patches as plt_patches
import numpy as np
# from scipy.optimize import curve_fit

import fit_functions as ff
import epid_peak_functions_rjnamm as epif
import datetime_functions as dtf

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "11.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"

YEAR_START = 1993
YEAR_END = 1993
MONTH_START, DAY_START = 7, 1  # July 1st
MONTH_END, DAY_END = 6, 30  # June 30th
CITIES = [
    {'name': 'Saint Petersburg', 'code': 'spb', 'marker': 'o'},
    {'name': 'Moscow', 'code': 'msk', 'marker': '*'},
    # {'name': 'Novosibirsk', 'code': 'nsk', 'marker': '^'},
]
# available matplotlib markers: "8><^vodDhH*ps'


def save_plot(city_name, data_list, data_list_biased, a1, a2, epidemic_curves,
              transit_curve1, transit_curve2):
    """Build and save one plot for a given dataset
    @param city_name: str, City's name
    @param data_list: np.ndarray, Dataset
    @param data_list_biased: np.ndarray, Dataset of biased values
    @param a1: np.ndarray with single lower non-flu ARI level
    @param a2: np.ndarray with single upper non-flu ARI level
    @param epidemic_curves: List[np.ndarray]
    @param transit_curve1: np.ndarray
    @param transit_curve2: np.ndarray
    """
    mfig = plt.figure(figsize=(17, 6.0))
    mfig.subplots_adjust(hspace=.3)

    fig1 = plt.subplot(111)
    fig1.set_title('ARI incidence, {0}, {1} to {2}'.format(
        city_name,
        dtf.convertfDateToFancyString(data_list[0, 0]),
        dtf.convertfDateToFancyString(data_list[-1, 0])),
        fontsize=24
    )

    days_list = data_list[..., 0]
    x_column = range(0, len(data_list), 1)
    x_labels_new = [dtf.returnProperDayNames(i) for i in days_list]
    plt.xticks(range(1, len(x_column) + 1, 1), x_labels_new)

    # create an array to mark actual data (a thursday for every week)
    x_thursday = dtf.returnThursdayMarks(data_list)

    incid_list_biased = data_list_biased[..., 3]
    plt.plot(x_column, incid_list_biased, 'c', label='Under-reported data', linewidth=2)
    plt.plot(x_thursday, incid_list_biased[x_thursday], 'co')

    incid_list = data_list[..., 3]
    epid_list = data_list[..., 4]
    plt.plot(x_column, incid_list, 'b', linewidth=2)
    plt.plot(x_thursday, incid_list[x_thursday], 'bo')

    # plotting epidemic curves
    epif.plotEpidCurves(epidemic_curves, days_list, x_thursday)
    epif.plotLevelTransitions(transit_curve1, transit_curve2, incid_list, x_thursday)

    x_data = np.linspace(0, len(incid_list) - 1, len(incid_list))
    plt.plot(x_data, ff.func(x_data, a1), "b--", label='Lower non-flu ARI level', linewidth=4)
    plt.plot(x_data, ff.func(x_data, a2), "r--", label='Higher non-flu ARI level', linewidth=4)
    plt.legend(loc='upper left', fontsize=18)

    plot_month_grid(fig1, data_list[0, 0], data_list[-1, 0])
    epif.plotEpidemicOutbreak(fig1, epid_list)

    # out_filename = '{0}-{1}.png'.format(year, year + 1)
    # plt.savefig(out_filename, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()


def plot_separate():
    """Plot epidemic curves of all cities for each year at the separate graphs"""
    for city in CITIES:
        for year in range(YEAR_START, YEAR_END + 1):
            date_start = datetime(year=year, month=MONTH_START, day=DAY_START).date()
            date_end = datetime(year=year + 1, month=MONTH_END, day=DAY_END).date()

            data_list_new = dtf.getDataArrayforGivenTimePeriod(
                np.loadtxt(r'out_dbases\\flu_dbase_%s.txt' % city['code']),
                date_start, date_end)
            data_list_biased = dtf.getDataArrayforGivenTimePeriod(
                np.loadtxt(r'out_dbases\\flu_dbase_%s_biased.txt' % city['code']),
                date_start, date_end)

            incidents = data_list_new[..., 3]
            epidemic_markers = data_list_new[..., 4]
            a1, a2, epidemic_peaks, epidemic_curves = epif.extractEpidStages(incidents, epidemic_markers)

            transit_curve1, transit_curve2 = epif.find_transit_curves(incidents, a1, a2)
            # ILI_high_beg_index = int(transit_curve1[-1]) if len(transit_curve1) > 0 else 0
            # ILI_high_end_index = int(transit_curve2[0]) if len(transit_curve2) > 0 else 0
            save_plot(city['name'], data_list_new, data_list_biased, a1, a2, epidemic_curves,
                      transit_curve1, transit_curve2)


def plot_epidemic_curves(city, days_list, x_thursday, epidemic_curves):
    """Plot epidemic curves for concrete city in time period. Assert matplotlib.plot is initialized
    @param city: dict(), 'name': Name of the city, 'marker': str() for matplotlib markers
    @param days_list: np.ndarray, Dataset
    @param x_thursday: List[int]
    @param epidemic_curves: List[np.ndarray]
    """
    for epidemic_curve in epidemic_curves:
        plt.plot(epidemic_curve[..., 0], epidemic_curve[..., 1], 'r', linewidth=2)
        days_list_peak = days_list[[int(i) for i in epidemic_curve[..., 0]]]

        epidemic_days = epidemic_curve[..., 0]

        x_thursday_peak = [int(i) for i in epidemic_days if i in x_thursday]
        x_thursday_peak_indices = [i for i in range(0, len(epidemic_days)) if epidemic_days[i] in x_thursday]

        plt.plot(x_thursday_peak, epidemic_curve[x_thursday_peak_indices, 1], 'r' + city['marker'], linewidth=2)
        filename_curve = 'epi_' + str(int(days_list[epidemic_curve[0, 0]])) + '.txt'
        np.savetxt(filename_curve, np.column_stack((days_list_peak, epidemic_curve[..., 1])), fmt="%d %d")


def plot_level_transitions(city, transit_curve1, transit_curve2, incidents, x_thursday):
    """Plot epidemic level transition curves for concrete city. Assert matplotlib.plot is initialized
    @param city: dict(), 'name': Name of the city, 'marker': str() for matplotlib markers
    @param transit_curve1: np.array
    @param transit_curve2: np.array
    @param incidents: np.ndarray
    @param x_thursday: List[int]
    """
    incidents = np.array(incidents)

    if len(transit_curve1):
        x_thursday_curve1 = [int(i) for i in transit_curve1 if i in x_thursday]
        plt.plot(transit_curve1, incidents[transit_curve1], "g", linewidth=2)
        plt.plot(x_thursday_curve1, incidents[x_thursday_curve1], "g" + city['marker'], linewidth=4)

    if len(transit_curve2):
        x_thursday_curve2 = [int(i) for i in transit_curve2 if i in x_thursday]
        plt.plot(transit_curve2, incidents[transit_curve2], "g", linewidth=2)
        plt.plot(x_thursday_curve2, incidents[x_thursday_curve2], "g" + city['marker'], linewidth=4)


def plot_month_grid(m_plot, date_start, date_finish):
    """Points out 1st day of each month in interval on incidence graph
    @type m_plot: matplotlib plot
    @param date_start: str, date formed by "%Y%m%d" string
    @param date_finish: str, date formed by "%Y%m%d" string
    """
    date_start = dtf.convertFloatToDate(date_start)
    date_finish = dtf.convertFloatToDate(date_finish)
    days_count = (date_finish-date_start).days

    for i in range(0, days_count+1, 1):
        if (date_start + timedelta(days=i)).day == 1:
            m_plot.axvline(i, color='gray', linestyle='dashed', linewidth=0.5)


def plot_at_same_graph():
    """Plot epidemic curves of all cities at the same graph for each year"""
    for year in range(YEAR_START, YEAR_END + 1):
        date_start = datetime(year=year, month=MONTH_START, day=DAY_START).date()
        date_end = datetime(year=year + 1, month=MONTH_END, day=DAY_END).date()

        mfig = plt.figure(figsize=(17, 6.0))
        mfig.subplots_adjust(hspace=.3)

        fig1 = plt.subplot(111)
        fig1.set_title('ARI incidence, {0} to {1}'.format(date_start, date_end), fontsize=24)
        legend_items = []

        for city in CITIES:
            data_list = dtf.getDataArrayforGivenTimePeriod(
                np.loadtxt(r'out_dbases\\flu_dbase_%s.txt' % city['code']),
                date_start, date_end)
            # data_list_biased = dtf.getDataArrayforGivenTimePeriod(
            #     np.loadtxt(r'out_dbases\\flu_dbase_%s_biased.txt' % city['code']),
            #     date_start, date_end)

            incidents = data_list[..., 3]
            epidemic_markers = data_list[..., 4]
            a1, a2, epidemic_peaks, epidemic_curves = epif.extractEpidStages(incidents, epidemic_markers)
            transit_curve1, transit_curve2 = epif.find_transit_curves(incidents, a1, a2)

            # create an array to mark actual data (a thursday for every week)
            x_thursday = dtf.returnThursdayMarks(data_list)

            days_list = data_list[..., 0]
            x_column = range(0, len(data_list), 1)
            x_labels_new = [dtf.returnProperDayNames(i) for i in days_list]
            plt.xticks(range(1, len(x_column) + 1, 1), x_labels_new)
            legend_items.append(plt_lines.Line2D([], [], color='black', markersize=10,
                                marker=city['marker'], label=city['name']))

            plt.plot(x_column, incidents, 'b', linewidth=2)
            plt.plot(x_thursday, incidents[x_thursday], 'b' + city['marker'])
            plot_level_transitions(city, transit_curve1, transit_curve2, incidents, x_thursday)
            plot_epidemic_curves(city, days_list, x_thursday, epidemic_curves)

        legend_items.extend([
            plt_patches.Patch(color='red', label='Epidemic outbreak'),
            plt_patches.Patch(color='green', label='Level transition'),
            plt_patches.Patch(color='blue', label='Reported data'),
        ])
        plt.legend(loc='upper left', numpoints=1, handles=legend_items, prop={'size': 16})
        plot_month_grid(fig1, data_list[0, 0], data_list[-1, 0])
        plt.show()
        plt.close()


if __name__ == '__main__':
    plot_separate()
    # plot_at_same_graph()
