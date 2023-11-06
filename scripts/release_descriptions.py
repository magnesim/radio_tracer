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

