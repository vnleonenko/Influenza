import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.interpolate import interp1d
import datetime
import calendar

#Converts weekly incidence data to daily and stores it in a new file

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
    return time_from_start/7 + 1 #week number

def fetchIncidenceList( filename ):
    #returns a list of weekly flu incidence with weeks and years
    #input file format: (year; week_num; inc_number; isEpidemic)

    reader = csv.reader(open(filename), delimiter=';')
    next(reader) #skipping the header
    my_list = list(reader)
    res_list=[row[0:4] for row in my_list] #taking only the necessary columns: (year; week_num; inc_number; epidemic_mark)
    return res_list

def extractDaysFromTimeline ( hourly_column ):
    #taking every fourth element of the list to extract dates
    daily_column = []
    for time_info in hourly_column[0::4]: #every 4th
        daily_column.append(str(int( ( int(time_info) / 100)) )) #removing hour time data

    daily_col_dates = [datetime.datetime.strptime(i, "%Y%m%d") for i in daily_column]

    return [i.date() for i in daily_col_dates]

def generateDates ( date_start, date_finish ):

    dates_arr = []

    #txtdate = '1.01.'+str(year)
    mydate = date_start #datetime.datetime.strptime(txtdate, "%d.%m.%Y")

    delta = datetime.timedelta(days=1)
    while mydate<date_finish:
        dates_arr.append(mydate)
        mydate+=delta

    M_res = np.column_stack((dates_arr))
    return M_res

def generateEpidMarkerArray( is_epid_nums, size ):
    res = np.zeros(size)

    is_epid_nums_corrected = [i for i in is_epid_nums if i<size and i>=0]
    res[is_epid_nums_corrected] = 1;
    return list(res)


def mergeIncidWithData ( res_list, city_mark ):
    #list structure: year; week; inc_num; epid_marker
    incidence_col_weeks = [int(row[2]) for row in res_list]

    year_first = int(res_list[0][0])
    week_first = int(res_list[0][1])
    year_last = int(res_list[len(res_list)-1][0])
    week_last = int(res_list[len(res_list)-1][1])

    date_start = fluinst_to_gregorian(year_first, week_first, 4)
    date_finish = fluinst_to_gregorian(year_last, week_last, 4)

    row_num = (date_finish-date_start).days+1 #from thursday to thursday

    print("Row_num:",row_num)

    x_interp = []
    is_epid_nums = [] #day numbers (from 0 to len) with epidemic markers

    for i in range(0,row_num,7):
        x_interp.append(i)
        if res_list[int(i/7)][3] == '1':
            for j in range(0,7,1):
                is_epid_nums.append(i-3+j) #from sunday to monday

    #generating the column of epidemic outbreak indicator
    is_epid = generateEpidMarkerArray( is_epid_nums, row_num )

    incidence_col_thursdays = [int(x/7.0) for x in incidence_col_weeks]

    print(incidence_col_weeks)
    print(incidence_col_thursdays)

    f = interp1d(x_interp,incidence_col_thursdays, kind='cubic')
    x_new = range(0,row_num, 1)
    inc_column_days = [int(i) for i in f(x_new)]


    #adding other data
    array_dates_part = generateDates(date_start, date_finish)
    print(array_dates_part)


    print('Now')
    print(len(array_dates_part))
    print(len(inc_column_days))
    print(len(is_epid))

    list_total = np.column_stack((array_dates_part,inc_column_days, is_epid))

    return list_total


years_space = range(1986,2016, 1)
city_mark_set = ['spb', 'msk', 'nsk']
prelim_treatment_type = ['orig', 'corrected']

for city_mark in city_mark_set:
    for treat in prelim_treatment_type:
        inc_list = fetchIncidenceList(r'input_raw_fixed\\zab_'+city_mark+'_' + treat + '.csv')

        res_list = mergeIncidWithData(inc_list, city_mark)
        fname = 'flu_inc_'+city_mark+'_'+treat+'.txt'
        np.savetxt(fname, np.column_stack(([datetime.datetime.strftime(i, "%Y%m%d") for i in res_list[...,0]],
                                                 res_list[...,1])), fmt="%s %d")


