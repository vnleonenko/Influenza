from datetime import datetime
from datetime import timedelta

def getDataArrayforGivenTimePeriod( data_array, date_start, date_finish ):
    "Extracting a part of data which corresponds to the given dates interval"

    data_col = [datetime.strptime(str(int(i)), '%Y%m%d').date() for i in data_array[...,0]]
    days = date_start-data_col[0]
    index_start = days.days
    days = data_col[-1] - date_finish
    index_finish=len(data_col)-days.days

    return data_array[index_start:index_finish]

def russianMonth(i):
    months = ['января', 'февраля', 'марта', 'апреля', 'мая', 'июня', 'июля', 'августа', 'сентября', 'октября', 'ноября', 'декабря']
    return months[i-1]

def convertfDateToFancyString( fdate ):
    mydate = datetime.strptime(str(int(fdate)), "%Y%m%d")
    return datetime.strftime(mydate, "%d %b %Y")

def convertfDateToFancyString_rus( fdate ):
    mydate = datetime.strptime(str(int(fdate)), "%Y%m%d")
    print(mydate)
    print(mydate.day)
    print(mydate.month)
    print(mydate.year)
    return str(mydate.day)+' '+russianMonth(mydate.month)+' '+str(mydate.year)+ ' г.' #datetime.strftime(mydate, "%d %b %Y")

def convertFloatToDate( fdate ):
    return datetime.strptime(str(int(fdate)), "%Y%m%d")

def returnProperDayNames( fdate ):
    #returns a day and month in a string format for every first day of the month
    mydate = datetime.strptime(str(int(fdate)), "%Y%m%d")
    if (mydate.day == 1):
        myperiod = datetime.strftime(mydate, "%d.%m")#"%d %b")
    else:
        myperiod = ' '
    return myperiod

def convertDateToStringDM( date ):
    return datetime.strftime(date, "%d %b")

def returnThursdayMarks( init_list ):
    "Returns marks on thursday days (ones taken from real data) for the array to be put in the plot"

    day_start = datetime.strptime(str(int(init_list[0][0])), "%Y%m%d")

    day_num = day_start.weekday() #days from 0 (Monday) to 6 (Sunday)
    if (day_num>3):
        thursday_index= 10-day_num
    else:
        thursday_index=3-day_num

    x_thursdays = []

    for i in range(thursday_index,len(init_list),7):
        x_thursdays.append(i)

    return x_thursdays

def plotSpecialDays ( plt, fdate_start, fdate_finish ):
    "Points out particular days on incidence graph"
    date_start = convertFloatToDate(fdate_start)
    date_finish = convertFloatToDate(fdate_finish)
    #storing holidays in int format mmdd
    holiday_set = ['09 Jan', '23 Feb', '08 Mar', '01 May', '02 May', '03 May', '04 May', '08 May', '09 May', '04 Nov', '12 Dec', '31 Dec', '01 Jan', '02 Jan', '03 Jan', '04 Jan', '05 Jan', '06 Jan', '07 Jan', '08 Jan', '10 Jan']

    ndays = (date_finish-date_start).days

    for i in range(0, ndays+1, 1):
        date_cur = date_start+timedelta(i)
        date_str = convertDateToStringDM(date_cur)
        #if date_str in holiday_set:
         #   fig1.axvline(i, color='green')

        if date_cur.day == 1:
            plt.axvline(i, color='gray', linestyle='dashed', linewidth=0.5)
            plt.axvline(i, color='gray', linestyle='dashed', linewidth=0.5)
            plt.axvline(i, color='gray', linestyle='dashed', linewidth=0.5)

