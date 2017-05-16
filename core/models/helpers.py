#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def max_elem_index(iterable):
    """returns the index of a highest incidence"""
    list_ = list(iterable)
    return list_.index(max(list_))


def calculate_tpeak_bias(y_model, data):
    peak_index_real = max_elem_index(data)
    peak_index_model = max_elem_index(y_model)
    delta = peak_index_model - peak_index_real
    return delta
