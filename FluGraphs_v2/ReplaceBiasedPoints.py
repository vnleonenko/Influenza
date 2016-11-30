import numpy as np
import csv
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import datetime

def fluinst_year_start(year):
    "The gregorian calendar date of the first day of the given year for the numeration of Flu Institute (1st Jan is always Week 1)"
    #Which means the first day of the first FluInst week
    first_jan = datetime.date(year, 1, 1) #is always contained in week 1
    delta = datetime.timedelta(first_jan.isoweekday()-1)
    return first_jan - delta

def fluinst_to_gregorian(year, week, day):
    "Gregorian calendar date for the given FluInst year, week and day"
    #days are from one to seven, thursday is fourth
    year_start = fluinst_year_start(year)
    return year_start + datetime.timedelta(days=day-1, weeks=week-1)

def gregorian_to_fluinst( date ):
    "Returns a week in FluInstutute numeration for the date given"
    year = date.year
    year_start = fluinst_year_start(year)
    time_from_start = date-year_start
    return time_from_start.days/7 + 1 #week number

#def week_to_fix ( date ):
#    "Uses gregorian_to_fluinst to find the week to fix for the holiday bias"
#    week_num = gregorian_to_fluinst(date)
#    if date.weekday()>3: #the holiday stands after thursday, which is current week incidence point
#        week_num+=1
#    return week_num

def findNPoints ( y_input, N ):
    "Find N points between four known ones using linear interpolation"
    x_input=[0,1,N+2,N+3] #0,1,2,3,4,5
    x = range(1,N+3)
    f=interp1d(x_input,y_input, kind='cubic')
    return f(x)

def replaceBiasPoints( inc_list_weeks ):
    "Replacing the biased points corresponding to holidays (2 points in Jan, one in Mar and Nov) by their linear interpolation"
    #Week nums: 1,2, one for March, 8th and one for Nov, 7th
    #inc_list_weeks contains year, week and inc_number
    #We suppose that the data is available for all the weeks 1..53 and sorted by week number

    week1_jan = 1
    #week2_jan = 2


    list_len = len([int(row[2]) for row in inc_list_weeks])


    for i in range(0, list_len):
        year = int(inc_list_weeks[i][0])
        week2_jan = gregorian_to_fluinst(datetime.date(year,12,31))
        week_feb = gregorian_to_fluinst(datetime.date(year,2,23))
        week_march = gregorian_to_fluinst(datetime.date(year,3,8))
        week1_may = gregorian_to_fluinst(datetime.date(year,5,1))
        week2_may = gregorian_to_fluinst(datetime.date(year,5,9))
        week_jun = gregorian_to_fluinst(datetime.date(year,6,12))
        week1_nov = gregorian_to_fluinst(datetime.date(year,11,4))
        week2_nov = gregorian_to_fluinst(datetime.date(year,11,7))
        week_dec = gregorian_to_fluinst(datetime.date(year,12,12))

        weeks = [week1_jan, week2_jan, week_feb, week_march, week1_may, week2_may, week1_nov, week2_nov, week_dec]

        cur_week = int(inc_list_weeks[i][1])

        if cur_week in weeks:
            if cur_week == week1_jan:
                incidence11 = int(inc_list_weeks[i-2][2])
                incidence12 = int(inc_list_weeks[i-1][2])
                incidence21 = int(inc_list_weeks[i+2][2])
                incidence22 = int(inc_list_weeks[i+3][2])

                y_out = findNPoints((incidence11, incidence12, incidence21, incidence22), 2)

                if (y_out[1]-float(inc_list_weeks[i][2])>0):
                    inc_list_weeks[i][2] = y_out[1]
                if (y_out[2]-float(inc_list_weeks[i+1][2])>0):
                    inc_list_weeks[i+1][2] = y_out[2]

            else:
                if cur_week == week2_jan:
                    incidence11 = int(inc_list_weeks[i-3][2])
                    incidence12 = int(inc_list_weeks[i-2][2])
                    incidence21 = int(inc_list_weeks[i+1][2])
                    incidence22 = int(inc_list_weeks[i+2][2])

                    y_out = findNPoints((incidence11, incidence12, incidence21, incidence22), 2)

                    if (y_out[1]-float(inc_list_weeks[i-1][2])>0): #correcting only lower values
                        inc_list_weeks[i-1][2] = y_out[1]

                    if (y_out[2]-float(inc_list_weeks[i][2])>0):
                        inc_list_weeks[i][2] = y_out[2]

                else:
                    incidence11 = int(inc_list_weeks[i-2][2])
                    incidence12 = int(inc_list_weeks[i-1][2])
                    incidence21 = int(inc_list_weeks[i+1][2])
                    incidence22 = int(inc_list_weeks[i+2][2])

                    y_out = findNPoints((incidence11, incidence12, incidence21, incidence22), 1)

                    if (y_out[1]-float(inc_list_weeks[i][2])>0):
                        inc_list_weeks[i][2] = y_out[1]

    return inc_list_weeks

def fetchIncidenceList( filename ):
    #returns a list of weekly flu incidence with weeks and years
    #input file format: (year; week_num; inc_number; isEpidemic)

    reader = csv.reader(open(filename), delimiter=';')
    next(reader) #skipping the head
    my_list = list(reader)
    res_list=[row[0:4] for row in my_list] #taking only the necessary columns: (year; week_num; inc_number; epidemic_mark)
    return res_list


def writeListToFile( res_list, filename):
    with open(filename, 'w') as file_:
        file_.write("Header \n")
        for item in res_list:
            #print(("%d; %d; %d; %s") % (int(item[0]), int(item[1]), int(item[2]), str(item[3])) )
            file_.write(("%d;%d;%d;%s\n") % (int(item[0]), int(item[1]), int(item[2]), str(item[3])) )

    return

list  = fetchIncidenceList('zab_spb_1986-2016.csv')
#list  = fetchIncidenceList('zab_spb_2011-2012.csv')
x = range(0,len(list))

#plt.plot(x,[int(row[2]) for row in list],'b')
list = replaceBiasPoints(list)
#plt.plot(x,[int(row[2]) for row in list],'r')
#plt.show()
writeListToFile(list,'zab_corrected_v3.csv')


