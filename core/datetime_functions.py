#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Auxiliary procedures for datetime operations
    based on BaroyanRvachev/datetime_functions.py
"""

from datetime import datetime


def convert_float_to_date(float_date):
    return datetime.strptime(str(int(float_date)), "%Y%m%d").date()


def returnSomeDaysNameFromDate(date):
    # print(type(date))
    # print(date)
    if date.day == 1 or date.day == 5 or date.day == 10 \
            or date.day == 15 or date.day == 20 or date.day == 25:
        my_period = datetime.strftime(date, "%d %b")
    else:
        my_period = ' '
    return my_period


def convertDateToStringDM(date):
    return datetime.strftime(date, "%d %b")


def convertDateToStringMY(date):
    return datetime.strftime(date, "%b %Y")
