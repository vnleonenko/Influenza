#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Converts weekly incidence data to daily and stores it in a new file"""

import csv
import datetime

import matplotlib.pyplot as plt
import matplotlib.lines as plt_lines
import matplotlib.patches as plt_patches
import numpy as np
from scipy.interpolate import interp1d

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "3.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"

# YEAR_START = 1985
# YEAR_END = 2015


def flu_institute_year_start(year):
    """The gregorian calendar date of the given year's first day
    for the numeration of Flu Institute (1st Jan is always Week 1)
    Which means the first day of the first FluInst week
    """
    first_jan = datetime.date(year, 1, 1)  # is always contained in week 1
    delta = datetime.timedelta(first_jan.isoweekday()-1)
    return first_jan - delta


def flu_institute_to_gregorian(year, week, day):
    """Gregorian calendar date for the given FluInst year, week and day
    days are from one to seven, thursday is fourth
    """
    year_start = flu_institute_year_start(year)
    return year_start + datetime.timedelta(days=day-1, weeks=week-1)


def gregorian_to_flu_institute(date):
    """Returns a week in FluInstutute numeration for the given date"""
    year = date.year
    year_start = flu_institute_year_start(year)
    time_from_start = date-year_start
    return time_from_start/7 + 1  # week number


def fetch_incidence_list(filename):
    """input file format: (year; week_num; inc_number; isEpidemic)
    returns a list of weekly flu incidence with weeks and years
    """
    with open(filename) as f:
        reader = csv.reader(f, delimiter=';')
        next(reader)  # skipping the header
        # taking only the necessary columns: (year; week_num; inc_number; epidemic_mark)
        res_list = [row[0:4] for row in list(reader)]
    return res_list


# def extractDaysFromTimeline ( hourly_column ):
#     #taking every fourth element of the list to extract dates
#     daily_column = []
#     for time_info in hourly_column[0::4]: #every 4th
#         daily_column.append(str(int( ( int(time_info) / 100)) )) #removing hour time data
#
#     daily_col_dates = [datetime.datetime.strptime(i, "%Y%m%d") for i in daily_column]
#
#     return [i.date() for i in daily_col_dates]


def generate_days_interval(date_start, date_finish):
    """Generate array of days from interval [date_start; date_end]
    @param date_start: datetime.datetime
    @param date_finish: datetime.datetime
    @return: np.array of datetime.datetime
    """
    dates_arr = []
    current_date = date_start

    while current_date <= date_finish:
        dates_arr.append(current_date)
        current_date += datetime.timedelta(days=1)

    return np.array(dates_arr)


def generate_epidemic_marker_list(epidemic_days_indexes, size):
    res = np.zeros(size)
    epidemic_days_corrected = [i for i in epidemic_days_indexes if 0 <= i < size]
    res[epidemic_days_corrected] = 1
    return list(res)


def merge_incidence_and_data(res_list):
    """list structure: year; week; inc_num; epidemic_marker"""
    incidence_col_weeks = [int(row[2]) for row in res_list]

    year_first = int(res_list[0][0])
    week_first = int(res_list[0][1])
    year_last = int(res_list[len(res_list)-1][0])
    week_last = int(res_list[len(res_list)-1][1])

    date_start = flu_institute_to_gregorian(year_first, week_first, 4)
    date_finish = flu_institute_to_gregorian(year_last, week_last, 4)

    row_num = (date_finish - date_start).days + 1  # from thursday to thursday
    print('Row_num: ', row_num)

    x_interp = []
    epidemic_day_indexes = []  # day numbers (from 0 to len) with epidemic markers

    for i in range(0, row_num, 7):
        x_interp.append(i)
        if res_list[int(i / 7)][3] == '1':
            for j in range(0, 7, 1):
                epidemic_day_indexes.append(i-3+j)  # from sunday to monday

    # generating the column of epidemic outbreak indicator
    is_epidemic = generate_epidemic_marker_list(epidemic_day_indexes, row_num)

    incidence_col_thursdays = [x // 7 for x in incidence_col_weeks]
    print('Incidence col weeks: ', incidence_col_weeks)
    print('Incidence col thursdays: ', incidence_col_thursdays)

    f = interp1d(x_interp, incidence_col_thursdays, kind='cubic')
    inc_column_days = [int(i) for i in f(range(0, row_num))]

    # adding other data
    days_array = generate_days_interval(date_start, date_finish)
    print("Days array", days_array)

    print('Now')
    print('len(days_array) = ', len(days_array))
    print('len(inc_column_days) = ', len(inc_column_days))
    print('len(is_epidemic) = ', len(is_epidemic))
    return np.column_stack((days_array, inc_column_days, is_epidemic, ))


def plot_week_grid(m_plot, date_start, date_finish):
    """Points out 1st day of each month in interval on incidence graph
    @type m_plot: matplotlib plot
    @param date_start: tuple, formed by ("%Y", "%W", ), i.e. year & week number
    @param date_finish: tuple, formed by ("%Y", "%W", ), i.e. year & week number
    """
    date_start = datetime.datetime.strptime(date_start[0] + '-W' + date_start[1] + '-4', "%Y-W%W-%w")
    date_finish = datetime.datetime.strptime(date_finish[0] + '-W' + date_finish[1] + '-4', "%Y-W%W-%w")
    weeks_count = (date_finish - date_start).days // 7

    for i in range(1, weeks_count + 1, 12):
        m_plot.axvline(i, color='gray', linestyle='dashed', linewidth=0.5)


def compare_week_data_with_interpolation():
    city_mark_set = ['spb', 'nsk']   # '['spb', 'msk', 'nsk']
    prelim_treatment_type = ['orig']  # , 'corrected']

    for city_mark in city_mark_set:
        for treat in prelim_treatment_type:
            inc_list = fetch_incidence_list(r'input_raw_fixed\\zab_%s_%s.csv' % (city_mark, treat))
            real_weeks = [int(row[2]) for row in inc_list]

            res_list = merge_incidence_and_data(inc_list)

            interpolated_weeks = [0]  # list of incident sums for consequent weeks
            current_week = res_list[0, 0].isocalendar()[1]
            for day, incidents, _ in res_list:
                week = day.isocalendar()[1]
                if week == current_week:
                    interpolated_weeks[-1] += incidents
                else:
                    interpolated_weeks.append(incidents)
                    current_week = week

            # interpolated_days = res_list[:-1, 1]
            # interpolated_weeks = np.sum(interpolated_days.reshape(-1, 7), axis=1)

            mfig = plt.figure(figsize=(17, 6.0))
            mfig.subplots_adjust(hspace=.3)
            fig1 = plt.subplot(111)
            fig1.set_title('Interpolation error for {0}, {1}'.format(
                city_mark, treat), fontsize=24)
            x_column = range(0, len(real_weeks))
            weeks = ['{} {}'.format(row[0], row[1]) for row in inc_list]

            # Cut the first and the last weeks
            x_column, weeks, real_weeks, interpolated_weeks = \
                x_column[1:-1], weeks[1:-1], real_weeks[1:-1], interpolated_weeks[1:-1]

            # Calculate the standard deviation
            deviation = sum(abs(r - i) ** 2 for r, i in zip(real_weeks, interpolated_weeks))
            print((deviation / len(interpolated_weeks)) ** (1/2))

            plt.xticks(x_column, weeks, rotation=60)
            plt.plot(x_column, real_weeks, 'r', label='Real weekly data', linewidth=2)
            plt.plot(x_column, interpolated_weeks, 'b', label='Interpolated data', linewidth=2)
            plt.legend(loc='upper left', fontsize=18)
            start = (inc_list[0][0], str(int(inc_list[0][1]) + 1), )
            finish = (inc_list[-1][0], str(int(inc_list[-1][1]) - 1), )
            plot_week_grid(fig1, start, finish)
            plt.show()
            plt.close()


def main():
    city_mark_set = ['nsk']   # '['spb', 'msk', 'nsk']
    prelim_treatment_type = ['orig', 'corrected']

    for city_mark in city_mark_set:
        for treat in prelim_treatment_type:
            inc_list = fetch_incidence_list(r'input_raw_fixed\\zab_%s_%s.csv' % (city_mark, treat))

            res_list = merge_incidence_and_data(inc_list)
            filename = 'flu_inc_%s_%s.txt' % (city_mark, treat)
            date_strings = ([datetime.datetime.strftime(i, "%Y%m%d") for i in res_list[..., 0]], res_list[..., 1])
            np.savetxt(filename, np.column_stack(date_strings), fmt="%s %d")

if __name__ == '__main__':
    compare_week_data_with_interpolation()
    # main()
