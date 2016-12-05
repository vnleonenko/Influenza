#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FluOptimizer implementations
    Group 1. numpy
        - 'Nelder-Mead'
        - 'Powell'
        - 'CG'
        - 'BFGS'
        - 'Newton-CG'
        - 'L-BFGS-B'
        - 'TNC'
        - 'COBYLA'
        - 'SLSQP'
    Group 2. Evolution programming
        - Genetic
"""

from functools import partial

from deap import base
from deap import creator
from deap import tools
import numpy as np
from scipy.optimize import minimize

from optimizer import FluOptimizer

__author__ = "Vasily Leonenko (vnleonenko@yandex.ru)"
__copyright__ = "Copyright 2016, ITMO University"
__version__ = "5.0"
__maintainer__ = "Nikita Seleznev (ne.seleznev@gmail.com)"


class BFGSOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='BFGS', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class LBFGSBOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='L-BFGS-B', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class SLSQPOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='SLSQP', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class NelderMeadOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, *args, **kwargs):
        result = minimize(function, minimize_params, method='Nelder-Mead')
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class PowellOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='Powell', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class CGOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='CG', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class TNCOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='TNC', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


class COBYLAOptimizer(FluOptimizer):
    def optimize(self, function, minimize_params, minimize_params_range):
        result = minimize(function, minimize_params, method='COBYLA', bounds=minimize_params_range)
        return result.fun, tuple(result.x)  # fit value, final bunch of optimal values


# noinspection PyPep8Naming
class GeneticOptimizer(FluOptimizer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        assert hasattr(self.params, 'POPULATION_SIZE')
        assert hasattr(self.params, 'CX_PROBABILITY')
        assert hasattr(self.params, 'MUT_PROBABILITY')
        assert hasattr(self.params, 'GENERATIONS_COUNT')

        # Types
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMin)

    # Operators
    def mutGaussian(self, individual, mu, sigma, indpb):
        K, I0 = individual
        if np.random.random() < indpb:
            K += np.random.normal(mu, sigma)
            K = max(min(K, self.params.K_RANGE[1]), self.params.K_RANGE[0])
        if np.random.random() < indpb:
            I0 += np.random.normal(mu, sigma)
            I0 = max(min(I0, self.params.I0_RANGE[1]), self.params.I0_RANGE[0])

        return [K, I0],

    def optimize(self, function, *args, **kwargs):

        # Initialization
        def evaluate(x):
            return function(x),

        K_min, K_max = self.params.K_RANGE
        I0_min, I0_max = self.params.I0_RANGE
        np.random.seed(42)
        init_params = zip(
            np.random.uniform(K_min, K_max, self.params.POPULATION_SIZE),
            np.random.uniform(I0_min, I0_max, self.params.POPULATION_SIZE)
        )
        generate_K_I0 = partial(next, init_params)

        toolbox = base.Toolbox()
        toolbox.register("individual", tools.initIterate, creator.Individual, generate_K_I0)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate", self.mutGaussian, mu=0,
                         sigma=(self.params.K_RANGE[1] - self.params.K_RANGE[0])/2,
                         indpb=0.1)
        toolbox.register("select", tools.selTournament, tournsize=3)
        toolbox.register("evaluate", evaluate)

        population = toolbox.population(n=self.params.POPULATION_SIZE)

        # Evaluate the entire population
        fitnesses = map(toolbox.evaluate, population)
        for ind, fitness in zip(population, fitnesses):
            ind.fitness.values = fitness

        for g in range(self.params.GENERATIONS_COUNT):
            # Select the next generation individuals
            offspring = toolbox.select(population, len(population))
            # Clone the selected individuals
            offspring = list(map(toolbox.clone, offspring))

            # Apply crossover and mutation on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if np.random.random() < self.params.CX_PROBABILITY:
                    toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            for mutant in offspring:
                if np.random.random() < self.params.MUT_PROBABILITY:
                    toolbox.mutate(mutant)
                    del mutant.fitness.values

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = map(toolbox.evaluate, invalid_ind)
            for ind, fitness in zip(invalid_ind, fitnesses):
                ind.fitness.values = fitness

            # The population is entirely replaced by the offspring
            population[:] = offspring

        fit_min, fit_arg = 10e100, (0, 0)
        fitnesses = map(toolbox.evaluate, population)
        for ind, fitness in zip(population, fitnesses):
            print(fitness[0], ind)
            if fitness[0] < fit_min:
                fit_min = fitness[0]
                fit_arg = ind
        print(fit_min, fit_arg)
        return fit_min, fit_arg
