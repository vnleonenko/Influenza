# Some old code

# def m(day):
#    """returns the  modifier for daily attendance of sick persons to healthcare units"""
#    switcher = {
#        0: 1.1,
#        1: 1.07,
#        2: 1.0,
#        3: 0.95,
#        4: 0.9,
#        5: 0.79,
#        6: 0.3,
#    }
#    return switcher.get(day, 0)  # zero infectivity by default

# def refine_data_from_raw(y):
#     """refines daily data using the statistics on daily attendance coefficients"""
#     y_refined = [y[i] * m(i % 7) for i in range(0,len(y))]
#
#     #replacing missed data by -1
#     arr = np.array(y_refined)
#     arr[arr<0]=-1
#     return list(arr)

# def convert_data_to_raw (y):
#     """converts the model daily data back to observed data"""
#     return [y[i] / m(i % 7) for i in range(0,len(y))]


# from itertools import cycle
# from matplotlib.font_manager import FontProperties
# lines = ["-", "--", "-.", ":"]
# line_cycler = cycle(lines)
# font = {'family': 'Verdana', 'weight': 'normal'}
# matplotlib.rc('font', **font)
#
# # ________ визуализация ________
# def DrawNewInfCases(DAYINF_TUP, WEEKINF_TUP):
#     DAYINF,DAYX = DAYINF_TUP
#     WEEKINF,WEEKX = WEEKINF_TUP
#
#     # ------- СТАТИСТИКА ПО ЗАБОЛЕВАНИЮ: -------
#     # ТЕКУЩЕЕ:
#     # S - общее число здоровых
#     # E - текущее число инфицированных (латентная стадия)
#     # I - текущее число инфицированных (активная стадия)
#     # T - общее число выздоровевших
#     # R - общее число умерших
#     f1 = plt.figure()
#     fontP = FontProperties()
#     fontP.set_size('small')
#
#     plt.plot(DAYX, DAYINF, next(line_cycler), label='Число заразившихся (день)', linewidth=2.0)
#     #plt.plot(WEEKX, WEEKINF, next(linecycler), label='Число заразившихся (неделя)', linewidth=2.0)
#
#     plt.grid()
#     plt.legend(loc='best',fancybox=True, shadow=True)
#     plt.xlabel('время')
#     plt.ylabel('популяция')
#     plt.title('Число заразившихся (день/неделя)')
#     plt.show()
