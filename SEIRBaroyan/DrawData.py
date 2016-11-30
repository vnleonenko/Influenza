# -*- coding: cp1251 -*-

import pylab as p
#from itertools import accumulate
from matplotlib.font_manager import FontProperties
from matplotlib import rc
from itertools import cycle

lines = ["-","--","-.",":"]
linecycler = cycle(lines)
font = {'family': 'Verdana','weight': 'normal'}
rc('font', **font)

# ________ визуализация ________
def DrawNewInfCases(DAYINF_TUP,WEEKINF_TUP):
    DAYINF,DAYX = DAYINF_TUP
    WEEKINF,WEEKX = WEEKINF_TUP

    # ------- СТАТИСТИКА ПО ЗАБОЛЕВАНИЮ: -------
    # ТЕКУЩЕЕ:
    # S - общее число здоровых
    # E - текущее число инфицированных (латентная стадия)
    # I - текущее число инфицированных (активная стадия)
    # T - общее число выздоровевших
    # R - общее число умерших    
    f1 = p.figure()
    fontP = FontProperties()
    fontP.set_size('small')
    
    p.plot(DAYX, DAYINF, next(linecycler), label='Число заразившихся (день)', linewidth=2.0)
    #p.plot(WEEKX, WEEKINF, next(linecycler), label='Число заразившихся (неделя)', linewidth=2.0)

    p.grid()
    p.legend(loc='best',fancybox=True, shadow=True)
    p.xlabel('время')
    p.ylabel('популяция')
    p.title('Число заразившихся (день/неделя)')
    p.show()