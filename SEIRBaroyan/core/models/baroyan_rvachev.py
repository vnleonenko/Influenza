#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SI-model of Baroyan-Rvachev origin with discrete time

Version history:
    * v2 - using numpy.optimize instead of iterating through k
    * v3 - including I0 into fitted params
    * v4 - considering bias in t_peak (iterating through possible value range)
    * v5 - refactor

    TODO Usage FluOptimizer
"""

import numpy as np

from .helpers import calculate_tpeak_bias

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "5.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"


class BaroyanParams:
    N = 3000  # epid duration
    T = 8  # disease duration for a single person

    SIZE = 25
    DISABLE_RANDOM = False
    RANDOM_SEED = 42  # Used if DISABLE_RANDOM is True
    K_RANGE = (1.02, 1.6)
    I0_RANGE = (1.0, 1.0)  # (0.1, 100)
    TPEAK_BIAS_RANGE = range(-7, 7)  # (-3, 3)

    @staticmethod
    def q(day):
        """returns infectivity by the day of infection - function q in B-R"""
        switcher = {
            2: 1,
            3: 0.9,
            4: 0.55,
            5: 0.3,
            6: 0.15,
            7: 0.05,
        }
        return switcher.get(day, 0)  # zero infectivity by default

    # FIXME Params for BaroyanGeneticOptimizer
    POPULATION_SIZE = 500
    CX_PROBABILITY = 0  # 0.5
    MUT_PROBABILITY = 0.2  # 0.2
    GENERATIONS_COUNT = 50


# noinspection PyPep8Naming
class FitFunction:
    def __init__(self, params, data, rho, tpeak_bias):
        self.params = params
        self.data = data
        self.rho = rho
        self.tpeak_bias = tpeak_bias

    @staticmethod
    def calculate_peak_bias(x, y):
        x_peak = max(x)
        y_peak = max(y)
        return abs(x_peak - y_peak)

    def sum_ill(self, y, t):
        """summing the cumulative infectivity of the infected on the moment t"""
        sum_ = 0

        for epid_day in range(0, self.params.T):
            if t - epid_day < 0:
                y_cur = 0
            else:
                y_cur = y[t - epid_day]

            # sum_ = sum_ + y[t-epid_day]*q(epid_day)
            sum_ += y_cur * self.params.q(epid_day)
        return sum_

    def make_simulation(self, alpha, lam, rho, I0):
        y = np.zeros((self.params.N + 1))
        x = np.zeros((self.params.N + 1))

        # initial data
        x[0] = alpha * rho
        y[0] = I0

        for t in range(0, self.params.N):
            y[t + 1] = lam * x[t] * self.sum_ill(y, t) / rho
            # print(y[t+1])
            x[t + 1] = x[t] - y[t + 1]
        return y

    @staticmethod
    def calculate_r(x, y, delta):
        """calculating the fitting coefficient r
        x is real data, y is modeled curve
        delta is the difference between the epidemic starts in real data and modeled curve
        """
        sum1 = 0
        sum2 = 0
        for i in range(delta, delta + len(x)):
            if x[i - delta] > 0 and y[i] > 0:  # do not consider absent data which is marked by -1
                sum1 += x[i - delta] * y[i]
                sum2 += pow(y[i], 2)
        return float(sum1) / float(sum2)

    @staticmethod
    def calculate_dist_squared(x, y, delta):
        """calculating the fitting coefficient r
        x is real data, y is modeled curve
        delta is the difference between the epidemic starts in real data and modeled curve
        """
        sum_ = 0
        for i in range(delta, delta + len(x)):
            if x[i - delta] > 0 and y[i] > 0:  # do not consider absent data which is marked by -1
                sum_ += pow(x[i - delta] - y[i], 2)
        return sum_

    def calculate_s(self, k):
        """calculating the parameter s to find the initial values of alpha and rho"""
        sum_ = 0
        T = self.params.T
        for tau in range(0, T):
            sum_ += pow(k, T - tau) * self.params.q(tau)

        return float(pow(k, T + 1)) / float(sum_)

    def find_model_fit(self, k, I0):
        """Launching the simulation for a given parameter value and aligning the result to model"""
        rho, tpeak_bias, data = self.rho, self.tpeak_bias, self.data
        s = self.calculate_s(k)
        alpha = 1
        lam = s
        # print(s)

        y_model = self.make_simulation(alpha, lam, rho, I0)

        # Aligning output by incidence peaks
        delta = calculate_tpeak_bias(y_model, data) + tpeak_bias  # adding possible peak moment bias

        #######################################################################
        if delta < 0:
            # sys.stderr.write("Model peak index is to the left of data peak!\n")
            return 10e10, [], -1, -1
        #######################################################################

        # Searching for the correction coefficient
        r = FitFunction.calculate_r(data, y_model, delta)
        # alpha = r
        # lam = s/r

        y_model = [r * y_item for y_item in y_model]  # adjusting model curve for the best fit
        # print(y_model)

        dist2 = FitFunction.calculate_dist_squared(data, y_model, delta)

        peak_bias = FitFunction.calculate_peak_bias(data, y_model)

        # len_data = len(FluOptimizer.data)
        # plt.plot(range(0, len_data), FluOptimizer.data)
        # plt.plot(range(0, len_data), y_model[:len_data])
        # plt.show()

        return dist2, y_model, delta, peak_bias

    def fit_function(self, params):
        k, I0 = params
        dist2, y_model, delta, peak_bias = self.find_model_fit(k, I0)
        return dist2


# noinspection PyPep8Naming
class AbstractBaroyanOptimizer:
    def __init__(self, data, population_quantity, params):
        self.data = data
        self.rho = population_quantity
        assert isinstance(params, BaroyanParams)
        self.params = params
        self.res2 = FitFunction.find_residuals(self.data)

        self.k_opt = None
        self.I0_opt = None
        self.R_square_opt = 0
        self.tpeak_bias_opt = None

    @staticmethod
    def find_residuals(data_list):
        """Finding the squares of residuals between the real data and it math expectation"""
        res = 0
        mean = np.mean(data_list)
        for item in data_list:
            res += pow(item - mean, 2)
        return res

    def optimize(self, function, minimize_params, minimize_params_range):
        """

        :param function:
        :param minimize_params: (K, I0)
        :param minimize_params_range:
        :return: (value, args) -- fit value, final bunch of optimal values
        """
        raise NotImplemented

    def _get_fit_function(self, tpeak_bias):
        return FitFunction(self.params, self.data, self.rho, tpeak_bias).fit_function

    def _get_model_fit(self):
        return FitFunction(self.params, self.data, self.rho, self.tpeak_bias_opt).find_model_fit

    def fit_one_outbreak(self):
        params_range = (self.params.K_RANGE, self.params.I0_RANGE)

        # k_opt, R_square_opt, I0_opt, tpeak_bias_opt = 0, 0, 0, 0

        # generating unifromly distributed init values for k
        K_min, K_max = self.params.K_RANGE
        I0_min, I0_max = self.params.I0_RANGE
        size = self.params.SIZE
        if self.params.DISABLE_RANDOM:
            np.random.seed(self.params.RANDOM_SEED)
        init_params = zip(
            np.random.uniform(K_min, K_max, size),
            np.random.uniform(I0_min, I0_max, size)
        )

        for params in init_params:
            for tpeak_bias_cur in self.params.TPEAK_BIAS_RANGE:
                function = self._get_fit_function(tpeak_bias_cur)

                value, args = self.optimize(function, params, params_range)
                k_cur, I0_cur = args

                R_square = 1 - value/self.res2

                # if peak_bias < peak_bias_opt:
                if R_square > self.R_square_opt or self.k_opt is None:
                    self.k_opt = k_cur
                    self.I0_opt = I0_cur
                    self.R_square_opt = R_square
                    self.tpeak_bias_opt = tpeak_bias_cur
                    # print(R_square_opt)
                    # print(peak_bias_opt)
            # print(k_opt, R_square)

    def get_results(self):
        find_model_fit = self._get_model_fit()

        dist2, y_model, delta, peak_bias = find_model_fit(self.k_opt, self.I0_opt)

        R_square_opt = 1 - dist2/self.res2

        return y_model, R_square_opt, self.k_opt, self.I0_opt, self.tpeak_bias_opt, delta
