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

