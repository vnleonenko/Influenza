#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Моделирование процесса распространения гриппа на основе системы ДУ
Идентификация параметров модели на основе данных измерений градиентным методом локального поиска

 -*- Flu Model -*-
    S - Susceptible (восприимчивые к болезни)
    E - Exposed (латентные - первые заразившиеся)
    I - Infected (инфицированные с признаками болезни)
    R - Recovered (выздоровевшие)
 -*- ФОРМА УРАВНЕНИЙ -*-
    dS/dt  = -k1*S*I              <-- текущее здоровых
    dE/dt =  +k1*S*I - k2*E       <-- текущее инфицированных латентных
    dI/dt  = +k2*E - k3*I         <-- текущее инфицированны
    dR/dt  = +k3*I                <--- ВЫЗДОРОВЕЛО ВСЕГО (ТОЛЬКО РОСТ):
 -*- СМЕРТНОСТЬ -*-
    dR/dt  = +k4*T                <--- УМЕРЛО ВСЕГО (ТОЛЬКО РОСТ) -- не используется

Version history:
    * v2 - corrected version (starting from the first measurement point if it's to the left of the model curve start
    * v3 - batch version to check the predictive force of the model
    * v4 - changed algorithm for fitting model curve to data - peak adjustment replaced by the discrete moving of the epid modeling starting moment
    * v5 - modified version to distinguish incomplete data without information about the peaks
    * v5b - starting level (coef = 0.9) is moved to fit function
    * v5cb - working on improved predictive ability
    * v6 - finding out the height and the date of an outbreak + fancy pics, rjnamm
    * v6b - still rjnamm, drawing S, E, I, R populations on the graph without baseline addition
    * v6_rjnamm - corrections according to the RJNAMM article review (new infected calc mechanism via E, R2 calculation routine)
    * v6_rjnamm_aN - experimental version using non-immune ratio instead of vert scaling and peak level to fit the curve
    * v8 - initial parameter variation added (like in _aN_intervals)
    * v9 refactor

"""

import numpy as np
import datetime
from core import datetime_functions as dtf

from core.models.helpers import calculate_tpeak_bias

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "5.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"


class SEIRParams:
    SIZE = 25
    DISABLE_RANDOM = False
    RANDOM_SEED = 42  # Used if DISABLE_RANDOM is True

    SHIFT_RANGE = (5, 54)  # unused
    K0_RANGE = (0.0001, 0.0001001)  # ???
    K1_RANGE = (0.0000001, 50.0)  # ???
    K2_RANGE = (0.39, 0.390001)  # ???
    K3_RANGE = (0.133, 0.1330001)  # ???
    S_RATIO_RANGE = (0.001, 1.0)  # ratio of susceptible from the total population
    SHIFT_COEFFICIENT_RANGE = (0.7, 1.0)


# noinspection PyPep8Naming
class SEIRFitting:
    def __init__(self, params, data, rho, shift):
        self.params = params
        self.data = data
        self.rho = rho
        self.shift = shift

    @staticmethod
    def get_infected(S, E):
        """Возвращает число зараженых за период (new cases), модель.
        S - результат моделирования."""
        INF = [-S[i + 1] + S[i] - E[i + 1] + E[i] for i in range(len(S) - 1)]
        return INF

    @staticmethod
    def make_simulation(initial_infected, k1, k2, k3, N):
        # ________ ДУ ________
        # N = N*3 # удвоим интервал моделирования ????

        t = np.linspace(0, N - 1, N)

        def dX_dt(X, t=0):  # Return parameters of the infected populations.
            S, E, I, T = X
            return np.array([-k1 * S * I,
                             +k1 * S * I - k2 * E,
                             +k2 * E - k3 * I,
                             +k3 * I])

        # ________ интегрируем ДУ ________
        from scipy import integrate
        # inInf = 0.01    # начальное число зараженных
        # переделать на zeros
        X0 = np.array([1. - initial_infected, initial_infected, 0.0, 0.0])  # initials conditions: S,E,I,T
        X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
        S, E, I, T = X.T

        return SEIRFitting.get_infected(S, E), S[1:], E[1:], I[1:], T[1:]

    @staticmethod
    def model_epidemic_duration(M_data):
        """finds the duration of epidemic from the data curve"""

        i = 0
        while M_data[i] * 20000 < 1 and i < len(M_data) - 1:
            i += 1

        epid_beg = i

        i = len(M_data) - 1

        while M_data[i] * 20000 < 1 and i > 0:
            i -= 1

        epid_end = i

        duration = max(0, epid_end - epid_beg)
        # print(duration)

        return duration

    @staticmethod
    def shift_model_to_data(DAY, M, shift, s_ratio, pop_total):
        """Отрезать лишние значения и произвести перенормировку.
            Shifting the model curve and cutting the points which are outside the data boundaries
        """

        # Calculating the prospected model epid length under this parameter set
        param = SEIRFitting.model_epidemic_duration(M)
        # print("Epi_duration: ", param)

        shiftM1 = shift  # столько отрезаем от начала модели (сдвиг модели влево)
        shiftM2 = len(DAY) + shift  # столько отрезаем от конца модели

        # print("Shifts: ", shiftM1, shiftM2, shiftM2 - shiftM1)

        M_new = M[shiftM1:shiftM2]
        # print(len(M))
        # print(len(M_new))

        # M_new = M_new * s_ratio * float(pop_total)
        M_new = [i * s_ratio * float(pop_total) for i in M_new]

        return DAY, M_new, param

    @staticmethod
    def shift_data_to_model(dates, DAY, M, shift, s_ratio, S, E, I, R, pop_total):
        # the same procedure modified for the sake of visualisation
        # (no model trajectories cut, hence the output arrays have different sizes)---

        # Moving data points to the right for the sake of plotting model vs data
        # Also scaling model incidence along with E and I
        N_zeros = len(M) - len(DAY)
        assert N_zeros >= 0  # considering the case when model data is always of greater length than the real data

        date_begin = dtf.convert_float_to_date(dates[0])
        date_end = dtf.convert_float_to_date(dates[-1])
        dates_new = []
        data_left = shift

        DAY_new = []

        DAY_new.extend([0] * shift)
        DAY_new.extend(DAY)
        data_right = len(DAY_new)
        DAY_new = list(DAY_new)

        DAY_new.extend([0] * (N_zeros - shift))

        date_begin_adj = date_begin - datetime.timedelta(days=int(shift))
        date_end_adj = date_end + datetime.timedelta(days=int(N_zeros - shift))
        M_new = M

        days_between1 = date_begin - date_begin_adj
        dates_range1 = [date_begin_adj + datetime.timedelta(days=x) for x in range(0, days_between1.days)]

        days_between2 = date_end_adj - date_end
        dates_range2 = [date_end + datetime.timedelta(days=x) for x in range(0, days_between2.days)]

        dates_new.extend(dates_range1)
        dates_new.extend([dtf.convert_float_to_date(x) for x in dates])
        dates_new.extend(dates_range2)

        # Calculating the old parameter scaling_coef for the comparison purposes
        maxMeas = max(DAY)  # максимальное значение (измерения)
        maxModel = max(M)  # максимальное значение (модель)
        scaling_coef = s_ratio * pop_total * maxModel / maxMeas
        # print(scaling_coef)

        # M_new = M_new * s_ratio * pop_total
        M_new = [i * s_ratio * float(pop_total) for i in M_new]

        # S_scaled = S* s_ratio * pop_total
        # E_scaled = E* s_ratio * pop_total
        # I_scaled = I* s_ratio * pop_total
        # R_scaled = R* s_ratio * pop_total

        # removing the redundant zeros to the right of M_new

        current_end = len(M_new) - 1

        min_cases_show = 10  # FIXME somehow

        while M_new[current_end] < min_cases_show and current_end > data_right:
            # the value of M_new is less than a fixed number
            current_end -= 1

        return (dates_new[:current_end + 1], DAY_new[:current_end + 1],
                M_new[:current_end + 1]), data_left, data_right, scaling_coef

    def fit_function(self, params):
        K0, K1, K2, K3, S_RATIO, SHIFT_COEFFICIENT = params

        # print(K)
        # процедура моделирования (incidence)
        y_model, S, E, I, R = SEIRFitting.make_simulation(K0, K1, K2, K3, len(self.data) * 3)

        shift = min(self.data) * SHIFT_COEFFICIENT

        data = self.data - shift

        DAY, M, modelEpidemicDuration = SEIRFitting.shift_model_to_data(data, y_model, self.shift, S_RATIO, self.rho)

        fit_sum_square_difference = sum(pow(DAY - M, 2))

        return fit_sum_square_difference


# noinspection PyPep8Naming
class AbstractSEIROptimizer:
    def __init__(self, initial_data, data, dates, population_quantity, params):
        self.initial_data = initial_data
        self.data = data
        self.dates = dates
        self.rho = population_quantity

        assert isinstance(params, SEIRParams)
        self.params = params
        self.res2 = self.find_residuals(self.data)

        self.k_opt = None  # tuple of optimal (K0, K1, K2, K3)
        self.s_ratio_opt = None
        self.shift_coefficient_opt = None
        self.shift_opt = None
        self.R_square_opt = 0

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

    def _get_fit_function(self, shift):
        return SEIRFitting(self.params, self.data, self.rho, shift).fit_function

    def fit_one_outbreak(self):
        params_range = (
            self.params.K0_RANGE, self.params.K1_RANGE, self.params.K2_RANGE, self.params.K3_RANGE,
            self.params.S_RATIO_RANGE, self.params.SHIFT_COEFFICIENT_RANGE)  # rjnamm_rev2 real_biol

        # generating unifromly distributed init values for k
        K1_min, K1_max = self.params.K1_RANGE
        S_RATIO_min, S_RATIO_max = self.params.S_RATIO_RANGE
        SHIFT_COEFFICIENT_min, SHIFT_COEFFICIENT_max = self.params.SHIFT_COEFFICIENT_RANGE
        size = self.params.SIZE

        if self.params.DISABLE_RANDOM:
            np.random.seed(self.params.RANDOM_SEED)
        init_params = zip(
            [0.0001] * size,  # Dummy K0
            np.random.uniform(K1_min, K1_max, size),
            [0.39] * size,  # Dummy K2
            [0.133] * size,  # Dummy K3
            np.random.uniform(S_RATIO_min, S_RATIO_max, size),
            np.random.uniform(SHIFT_COEFFICIENT_min, SHIFT_COEFFICIENT_max, size)
        )
        init_params = list(init_params)

        for idx, params in enumerate(init_params):
            # FIXME Vasily's hot heuristics
            if self.R_square_opt and self.R_square_opt > 0.8 and idx > 1:
                print('Break!')
                break

            for shift in range(self.params.SHIFT_RANGE[0], self.params.SHIFT_RANGE[1]):
                function = self._get_fit_function(shift)

                value, args = self.optimize(function, params, params_range)
                K0, K1, K2, K3, S_RATIO, SHIFT_COEFFICIENT = args

                R_square = 1 - value/self.res2

                # if peak_bias < peak_bias_opt:
                if R_square > self.R_square_opt or self.k_opt is None:
                    self.k_opt = (K0, K1, K2, K3)
                    self.s_ratio_opt = S_RATIO
                    self.shift_coefficient_opt = SHIFT_COEFFICIENT
                    self.shift_opt = shift
                    self.R_square_opt = R_square

                    # FIXME Vasily's hot heuristics
                    if R_square > 0.9 and idx > 1:
                        print('Break!')
                        break

    def get_results(self):

        def calculate_peak_bias(x, y):  # FIXME if used at least twice
            return max(y) / max(x)

        # launching simulation with optimal param set for visualization purpose
        k0, k1, k2, k3 = self.k_opt
        y_model, S, E, I, R = SEIRFitting.make_simulation(
            k0, k1, k2, k3, len(self.initial_data) * 3)

        level_zero = min(self.data)
        vertical_shift = level_zero * self.shift_coefficient_opt

        model_data, data_left, data_right, scaling_coefficient = SEIRFitting.shift_data_to_model(
            self.dates, self.data - vertical_shift, y_model, self.shift_opt, self.s_ratio_opt,
            S, E, I, R, self.rho)

        peak_bias = calculate_peak_bias(self.initial_data + vertical_shift, model_data[2] + vertical_shift)
        tpeak_bias = calculate_tpeak_bias(self.initial_data + vertical_shift, model_data[2] + vertical_shift)

        blob = {
            "y_model": y_model,
            "R_square_opt": self.R_square_opt,
            "k_opt": self.k_opt,
            "s_ratio_opt": self.s_ratio_opt,

            "shift_coefficient_opt": self.shift_coefficient_opt,
            "shift_opt": self.shift_opt,

            "S": S, "E": E, "I": I, "R": R,

            "peak_bias": peak_bias,
            "tpeak_bias": tpeak_bias,

            "model_data": model_data,
            "data_left": data_left,
            "data_right": data_right,
            "scaling_coefficient": scaling_coefficient
        }
        return blob
