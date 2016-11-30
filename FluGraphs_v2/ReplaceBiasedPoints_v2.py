#v1 (without version num) - holidays are put manually and are the same for all the years
#v2 (Yulia's modification) - holidays are uploaded from the list and are subjected to change over the years

import csv
from scipy.interpolate import interp1d
import datetime
import math

def fluinst_year_start(year):
    "The gregorian calendar date of the first day of the given year for the numeration of Flu Institute (1st Jan is always Week 1)"
    #Which means the first day of the first FluInst week
    first_jan = datetime.date(year, 1, 1) #is always contained in week 1
    delta = datetime.timedelta(first_jan.isoweekday()-1)
    return first_jan - delta

def gregorian_to_fluinst( date ):
    "Returns a week in FluInstutute numeration for the date given"
    year = date.year
    year_start = fluinst_year_start(year)
    time_from_start = date-year_start
    return math.trunc(time_from_start.days/7 + 1) #week number

def findNPoints ( y_input, N ):
    "Find N points between four known ones using linear interpolation"
    x_input=[0,1,N+2,N+3] #0,1,2,3,4,5
    x = range(1,N+3)
    f=interp1d(x_input,y_input, kind='cubic')
    return f(x)

def replaceBiasPoints( inc_list_weeks, holidays ):
    "Replacing the biased points corresponding to holidays by their linear interpolation"
    #inc_list_weeks contains year, week and inc_number
    #We suppose that the data is available for all the weeks 1..53 and sorted by week number
    #holidays contains year, month and day

    week1_jan = 1
    week2_jan = 2

    list_len = len([int(row[2]) for row in inc_list_weeks])

    for year in range(1985, 2016):
        weeks = []
        for date in holidays:
            if int(date[0]) == year:
                week = gregorian_to_fluinst(datetime.date(year, int(date[1]), int(date[2])))
                if week not in weeks:
                    weeks.append(week)

        for i in range(0, list_len):

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


def writeListToFile( res_list, filename ):
    with open(filename, 'w') as file_:
        file_.write("Header \n")
        for item in res_list:
            file_.write(("%d;%d;%d;%s\n") % (int(item[0]), int(item[1]), int(item[2]), str(item[3])) )
    return

def readFromCsvToList( filename ):
    #return list with all data from csv file
    reader = csv.reader(open(filename), delimiter=';')
    next(reader)
    res_list = list(reader)
    return res_list

holidays = readFromCsvToList('holidays_1986-2016.csv')

list = fetchIncidenceList('zab_spb_1986-2016.csv')

list = replaceBiasPoints(list, holidays)

writeListToFile(list,'zab_corrected.csv')