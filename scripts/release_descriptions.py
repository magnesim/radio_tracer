import numpy as np
import pandas as pd
from datetime import datetime, timedelta



def releases_sellafield(dates=None, isotop=None, number=0):
    if not dates==None:
        d0 = dates[0]
        d1 = dates[1]
    else:
        d0 = datetime(1990,1,1)
        d1 = datetime(2023,1,1)

    if number>0:
        date_arr = pd.date_range(start=d0, end=d1, periods=number)
    else:
        date_arr = np.arange(d0,d1,timedelta(days=1))


    nx = len(date_arr)

    print('number',number,'nx',nx)

    if isotop=='129I':
        lbd = 1.3
        alpha = 1.
        zeta = 200
    elif isotop=='236U':
        lbd = -3.
        alpha = 1300.
        zeta = 300

    rel1 = zeta + alpha* np.exp(lbd * np.arange(nx)/float(nx))

    return date_arr, rel1



def releases_lahague(dates=None, isotop=None, number=0):
    date_arr, rel1 = releases_sellafield(dates=dates, isotop=isotop, number=number)
    if isotop=='129I':
        rel1 = rel1*20
    elif isotop=='236U':
        rel1 = rel1 *.1

    return date_arr, rel1





def monthly_release(isotop=None, dates=None, location=None):
    import pandas as pd

    folder = '../data'

    if not isotop in ['129I', '236U']:
        print(' Illegal isotope: ',isotop)
        print('Legal isotopes: [129I, 236U]')
        exit()

    if not location in ['Sellafield', 'LaHague']:
        print(' Illegel location: ',location)
        print('Legal locations: [Sellafield, LaHague]')
        exit()


    fn = '{}/{}_{}.csv'.format(folder, isotop, location.strip())

    # Read input file
    print('read from file:',fn)
    try: 
        df = pd.read_csv(fn, header=0, usecols=["Dato (YYYYMMDD)", "Number of Atoms"])
    except Exception:
        print('File not available: ',fn )
        exit()
    date_arr = np.array([datetime.strptime(str(item), '%Y%m%d') for item in df["Dato (YYYYMMDD)"] ])
    rel1 = np.array([item for item in df["Number of Atoms"]])

    date_arr = date_arr[np.where((date_arr>=dates[0]) & (date_arr<=dates[1]))]    
    rel1 = rel1[np.where((date_arr>=dates[0]) & (date_arr<=dates[1]))]



    return date_arr, rel1 




def get_daily_weights(release, dates):
    '''
    release:     [datetime objects, release magnitude]
    dates:       array of dates to be evaluated
    ntraj:       number of trajectories (int)
    '''
    
    import calendar

    [release_t, release_y] = release

    natoms_total = np.sum(release_y)
    rel_release = release_y / natoms_total 

    #ndays = int((dates[1] - dates[0]).total_seconds() / 86400 +1)
    ndays = int(len(dates))
    print(f'ndays: {ndays}',ndays)


#    date_arr = [dates[0] + timedelta(days=item) for item in range(ndays)]
    date_arr = dates

    
    weightsDay=np.zeros(ndays)
    release_year = np.array([item.year for item in release_t])
    release_month = np.array([item.month for item in release_t])


    for ii in range(ndays):
        yy = date_arr[ii].year
        mm = date_arr[ii].month 
        dim = calendar.monthrange(yy, mm)[1]
        try:
            tt = np.where((release_month == mm) & (release_year==yy) )[0][0]
        except Exception:
            print('Error: get_daily_weights: ', yy, mm)
            exit()
        weightsDay[ii] = rel_release[tt] / dim 


    return date_arr, weightsDay, natoms_total




def get_seed_date(file):
    import opendrift 
    from netCDF4 import Dataset, num2date

    nc = Dataset(file, 'r')
    age = nc.variables['age_seconds'][:]
    time = nc.variables['time'][:]
    time_unit = nc.variables['time'].units 
    ntra = len(nc.variables['trajectory'][:])
    nc.close()


    idx = np.argmin(age[:], axis=1)

    mindate=np.zeros(ntra)
    for tr in range(ntra):
        mindate[tr] = time[idx[tr]]
        
    mindate = np.array([num2date(item, time_unit, only_use_cftime_datetimes=False) for item in mindate])

    timelims = [num2date(time[0],time_unit, only_use_cftime_datetimes=False),
                num2date(time[-1],time_unit, only_use_cftime_datetimes=False)]

    return mindate, timelims, ntra



'''
def mon_rel_LH(isotop=None, dates=None):
    import pandas as pd

    folder = '../data'

    if isotop=='129I':
        fn = folder+'/129I_LaHague.csv'
    elif isotop=='236U':
        fn = folder+'/Modellfil-236ULaHague.csv'
    elif isotop==None:
        exit()

    # Read input file
 
    print('read from file:',fn)
    df = pd.read_csv(fn, header=0, usecols=["Dato (YYYYMMDD)", "Number of Atoms"])
    date_arr = np.array([datetime.strptime(str(item), '%Y%m%d') for item in df["Dato (YYYYMMDD)"] ])
    rel1 = np.array([item for item in df["Number of Atoms"]])

    date_arr = date_arr[np.where((date_arr>=dates[0]) & (date_arr<=dates[1]))]    
    rel1 = rel1[np.where((date_arr>=dates[0]) & (date_arr<=dates[1]))]

    return date_arr, rel1 


def mon_rel_SF(isotop=None, dates=None):
    import pandas as pd

    folder = '../data'

    if isotop=='129I':
        fn = folder+'/129I_Sellafield.csv'
    elif isotop=='236U':
        fn = folder+'/Modellfil-236USellafield.csv'
    elif isotop==None:
        exit()

    # Read input file
 
    print('read from file:',fn)
    df = pd.read_csv(fn, header=0, usecols=["Dato (YYYYMMDD)", "Number of Atoms"])
    date_arr = np.array([datetime.strptime(str(item), '%Y%m%d') for item in df["Dato (YYYYMMDD)"] ])
    rel1 = np.array([item for item in df["Number of Atoms"]])

    date_arr = date_arr[np.where((date_arr>=dates[0]) & (date_arr<=dates[1]))]    
    rel1 = rel1[np.where((date_arr>=dates[0]) & (date_arr<=dates[1]))]


    return date_arr, rel1 

'''
