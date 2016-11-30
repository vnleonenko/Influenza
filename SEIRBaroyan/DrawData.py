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

# ________ ������������ ________
def DrawNewInfCases(DAYINF_TUP,WEEKINF_TUP):
    DAYINF,DAYX = DAYINF_TUP
    WEEKINF,WEEKX = WEEKINF_TUP

    # ------- ���������� �� �����������: -------
    # �������:
    # S - ����� ����� ��������
    # E - ������� ����� �������������� (��������� ������)
    # I - ������� ����� �������������� (�������� ������)
    # T - ����� ����� �������������
    # R - ����� ����� �������    
    f1 = p.figure()
    fontP = FontProperties()
    fontP.set_size('small')
    
    p.plot(DAYX, DAYINF, next(linecycler), label='����� ������������ (����)', linewidth=2.0)
    #p.plot(WEEKX, WEEKINF, next(linecycler), label='����� ������������ (������)', linewidth=2.0)

    p.grid()
    p.legend(loc='best',fancybox=True, shadow=True)
    p.xlabel('�����')
    p.ylabel('���������')
    p.title('����� ������������ (����/������)')
    p.show()