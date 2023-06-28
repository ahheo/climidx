#!/usr/bin/env python3

from iris.util import mask_cube

from .func_ import *

# INDEX DICT: FORMAT:
#     NAME: (#, tInt_i, VARIABLES, FUNCTION, tInt_o, METADATA, GROUPS)
# where tInt_i, tInt_o are input and output temporal resolution.
# acceptable values and desciption:
#     tInt_i: ('mon', 'day', '6hr', '3hr', '1hr')
#     VARIABLES: input variables (using the same names as CORDEX/CMIPX)
#     FUNCTION (most fr. climidx module):
#         None or string: using pSTAT_cube() fr. module uuuu;
#                         see _d0() in _exe() in exe_.py
#         function type1: with cube(s) as input(s); see _d1() in _exe() in 
#                         exe_.py
#         function type2: with 1d array (np) as input(s); see _d2() in _exe()
#                         in exe_.py
#     tInt_o: ('month', 'season', 'year', 'hour-month', 'hour-season', 'hour')
#         PS: 'month' is exactly 'month-year'; similar for 'season',
#             'hour-month', 'hour-season', and 'hour'
#     EXE: (
#         function type: 0, 1, 2, and 3, ..., 8 based on 0-2 but deeling with
#                        some complicate situations
#         input variables to FUNCTION: None: same as VARIABLES; otherwise
#         positional and keyword arguments
#          )
#     METADATA:
#         a diction:
#             name: try standard name first; if failed long_name instead
#             units:
#             attrU: some notes
#         None: nothing modified but following the input nc files
#     GROUPS (for batch application;
#             one or more of the following;
#             free to add more):
#         'p': precipitation related
#         't': temperature related
#         'w': wind related
#         'r': radiation related
#         'c': consecutive days
i__ = {
    'ET': (
        1,
        'mon',
        ['evspsbl'],
        None,
        ['year'],
        (0, None),
        None,
        'p'),
    'HumiWarmDays': (
        2,
        'day',
        ['hurs', 'tas'],
        dHumiWarmDays_cube,
        ['season'],
        (1, None, (), {}),
        dict(
            name='humid warm days',
            units='days',
            attrU={'CLIMI': 'days with hurs > 90% & tas > 10 degree C'}
            ),
        'ht'),
    'EffPR': (
        3,
        'mon',
        ['evspsbl', 'pr'],
        None,
        ['season', 'year'],
        (0, ['pMINUSe']),
        dict(
            name='effective precipitation',
            attrU={'CLIMI': 'pr - evspsbl'}
            ),
        'p'),
    'PR7Dmax': (
        4,
        'day',
        ['pr'],
        dPr7_,
        ['year'],
        (2, None, 'y', {}, {}),
        dict(
            name='max 7-day precipitation',
            units='days',
            attrU={'CLIMI': 'maximum of rolling 7-day mean precipitation'}),
        'p'),
    'LnstDryDays': (
        5,
        'day',
        ['pr'],
        dLongestDryDays_,
        ['season'],
        (2, None, 'y', {}, {}),
        dict(
            name='longest dry days',
            units='days',
            attrU={'CLIMI': 'longest dry spell defined as consecutive days '
                            'with daily pr < 1 mm'}
            ),
        'p'),
    'DryDays': (
        6,
        'day',
        ['pr'],
        dDryDays_,
        ['year'],
        (1, None, (), {}),
        dict(
            name='dry days',
            units='days',
            attrU={'CLIMI': 'total days with daily pr < 1 mm'}
            ),
        'p'),
    'PR': (
        7,
        'mon',
        ['pr'],
        None,
        ['season', 'year'],
        (0, None),
        None,
        'p'),
    'PRmax': (
        8,
        'day',
        ['pr'],
        'MAX',
        ['month', 'year'],
        (0, None),
        dict(name='maximum precipitation rate'),
        'p'),
    'PRRN': (
        9,
        'mon',
        ['pr', 'prsn'],
        None,
        ['season', 'year'],
        (0, ['prrn']),
        dict(name='rainfall flux'),
        'p'),
    'PRSN': (
        10,
        'mon',
        ['prsn'],
        None,
        ['season', 'year'],
        (0, None),
        None,
        'p'),
    'NetRO': (
        11,
        'mon',
        ['mrro'],
        None,
        ['year'],
        (0, ['mrro_amjjas']),
        dict(attrU={'CLIMI': 'season: amjjas'}),
        'p'),
    'SD': (
        12,
        'mon',
        ['sund'],
        None,
        ['year'],
        (0, None),
        None,
        'r'),
    'RLDS': (
        13,
        'mon',
        ['rlds'],
        None,
        ['season'],
        (0, None),
        None,
        'r'),
    'RSDS': (
        14,
        'mon',
        ['rsds'],
        None,
        ['season'],
        (0, None),
        None,
        'r'),
    'CoolingDegDay': (
        15,
        'day',
        ['tasmax'],
        dDegDay_,
        ['month', 'year'],
        (1, None, (), dict(thr=20)),
        dict(
            name='degree day for cooling',
            attrU={'CLIMI': 'degree day for tasmax > 20 degree C'}
            ),
        't'),
    'ConWarmDays': (
        16,
        'day',
        ['tasmax'],
        dConWarmDays_,
        ['year'],
        (2, None, 'y', {}, {}),
        dict(
            name='consecutive warm days',
            units='days',
            attrU={'CLIMI': 'longest warm spell defined as consecutive days '
                            'with tasmax > 20 degree C'}
            ),
        'tc'),
    'TX': (
        17,
        'mon',
        ['tasmax'],
        None,
        ['year'],
        (0, None),
        None,
        't'),
    'WarmDays': (
        18,
        'day',
        ['tasmax'],
        dWarmDays_,
        ['season', 'year'],
        (1, None, (), {}),
        dict(
            name='warm days',
            units='days',
            attrU={'CLIMI': 'total days with tasmax > 20 degree C'}
            ),
        't'),
    'ColdDays': (
        19,
        'day',
        ['tasmax'],
        dColdDays_,
        ['season', 'year'],
        (1, None, (), {}),
        dict(
            name='cold days',
            units='days',
            attrU={'CLIMI': 'total days with tasmax < -7 degree C'}
            ),
        't'),
    'DegDay20': (
        20,
        'day',
        ['tas'],
        dDegDay_,
        ['year'],
        (1, None, (), dict(thr=20)),
        dict(
            name='degree day warmer than 20 degree C',
            attrU={'CLIMI': 'degree day for tas > 20 degree C'}
            ),
        't'),
    'DegDay8': (
        21,
        'day',
        ['tas'],
        dDegDay8_vegSeason_,
        ['year'],
        (2, None, 'y', {}, {}),
        dict(
            name='degree day under vegetation season',
            attrU={'CLIMI': 'degree day for tas > 8 degree C under '
                            'vegetation season (tas > 5 degree C)'}
            ),
        't'),
    'DegDay17': (
        22,
        'day',
        ['tas'],
        dDegDay_,
        ['year'],
        (1, None, (), dict(thr=17, left=True)),
        dict(
            name='degree day cooler than 17 degree C',
            attrU={'CLIMI': 'degree day for tas < 17 degree C'}
            ),
        't'),
    'TN': (
        23,
        'mon',
        ['tasmin'],
        None,
        ['year'],
        (0, None),
        None,
        't'),
    'SpringFrostDayEnd': (
        24,
        'day',
        ['tasmin'],
        dEndSpringFrost_,
        ['year'],
        (2, None, 'yd', {}, {}),
        dict(
            name='end of spring frost',
            units=1,
            attrU={'CLIMI': 'last day with tasmin < 0 degree C; '
                            'no later than 213.'}
            ),
        't'),
    'FrostDays': (
        25,
        'day',
        ['tasmin'],
        dFrostDays_,
        ['season', 'year'],
        (1, None, (), {}),
        dict(
            name='frost days',
            units='days',
            attrU={'CLIMI': 'total days with tasmin < 0 degree C'}
            ),
        't'),
    'TropicNights': (
        26,
        'day',
        ['tasmin'],
        dTropicNights_,
        ['year'],
        (1, None, (), {}),
        dict(
            name='tropical nights',
            units='days',
            attrU={'CLIMI': 'total days with tasmin > 17 degree C'}
            ),
        't'),
    'ZeroCrossingDays': (
        27,
        'day',
        ['tasmax', 'tasmin'],
        dZeroCrossingDays_cube,
        ['season'],
        (1, None, (), {}),
        dict(
            name='zero-crossing days',
            units='days',
            attrU={'CLIMI': 'days with tasmin < 0 degree C and '
                            'tasmax > 0 degree C'}
            ),
        't'),
    'VegSeasonDayEnd-5': (
        28,
        'day',
        ['tas'],
        dStartEndVegSeason_,
        ['year'],
        (3, None, 'yd', dict(thr=5)),
        dict(
            name='vegetation season (5) end',
            units=1,
            attrU={'CLIMI': 'vegetation season: from first day to last day '
                            'with4-day rolling mean tas > 5 degree C'}
            ),
        't'),
    'VegSeasonDayEnd-2': (
        29,
        'day',
        ['tas'],
        dStartEndVegSeason_,
        ['year'],
        (3, None, 'yd', dict(thr=2)),
        dict(
            name='vegetation season (2) end',
            units=1,
            attrU={'CLIMI': 'vegetation season: from first day to last day '
                            'with4-day rolling mean tas > 2 degree C'}
            ),
        't'),
    'VegSeasonDayStart-5': (
        30,
        'day',
        ['tas'],
        dStartEndVegSeason_,
        ['year'],
        (3, None, 'yd', dict(thr=5)),
        dict(
            name='vegetation season (5) start',
            units=1,
            attrU={'CLIMI': 'vegetation season: from first day to last day '
                            'with4-day rolling mean tas > 5 degree C'}
            ),
        't'),
    'VegSeasonDayStart-2': (
        31,
        'day',
        ['tas'],
        dStartEndVegSeason_,
        ['year'],
        (3, None, 'yd', dict(thr=2)),
        dict(
            name='vegetation season (2) start',
            units=1,
            attrU={'CLIMI': 'vegetation season: from first day to last day '
                            'with4-day rolling mean tas > 2 degree C'}
            ),
        't'),
    'VegSeasonLength-5': (
        32,
        'day',
        ['tas'],
        dStartEndVegSeason_,
        ['year'],
        (3, None, 'yd', dict(thr=5)),
        dict(
            name='vegetation season (5) length',
            units=1,
            attrU={'CLIMI': 'vegetation season: from first day to last day '
                            'with4-day rolling mean tas > 5 degree C'}
            ),
        't'),
    'VegSeasonLength-2': (
        33,
        'day',
        ['tas'],
        dStartEndVegSeason_,
        ['year'],
        (3, None, 'yd', dict(thr=2)),
        dict(
            name='vegetation season (2) length',
            units=1,
            attrU={'CLIMI': 'vegetation season: from first day to last day '
                            'with4-day rolling mean tas > 2 degree C'}
            ),
        't'),
    'SfcWind': (
        34,
        'day',
        ['sfcWind', 'uas', 'vas'],
        None,
        ['month', 'season', 'year'],
        (0, ['sfcWind']),
        None,
        'w'),
    'WindGustMax': (
        35,
        'day',
        ['wsgsmax'],
        'MAX',
        ['year'],
        (0, None),
        dict(name='maximum wind gust'),
        'w'),
    'WindyDays': (
        36,
        'day',
        ['wsgsmax'],
        dWindyDays_,
        ['season', 'year'],
        (1, None),
        dict(
            name='windy days',
            units='days',
            attrU={'CLIMI': 'total days with wsfsmax > 21 m s-1'}
            ),
        'w'),
    'PRgt10Days': (
        37,
        'day',
        ['pr'],
        dExtrPrDays_,
        ['season', 'year'],
        (1, None, (), {}),
        dict(
            name='days with heavy precipitation',
            units='days',
            attrU={'CLIMI': 'total days with daily pr > 10 mm'}
            ),
        'p'),
    'PRgt25Days': (
        38,
        'day',
        ['pr'],
        dExtrPrDays_,
        ['season', 'year'],
        (1, None, (), dict(thr=25)),
        dict(
            name='days with extreme precipitation',
            units='days',
            attrU={'CLIMI': 'total days with daily pr > 25 mm'}
            ),
        'p'),
    'SncDays': (
        39,
        'day',
        ['snc'],
        dSncDays_,
        ['year'],
        (1, None, (), {}),
        dict(
            name='days with snow cover',
            units='days',
            attrU={'CLIMI': 'total days with snc > 0'}
            ),
        'p'),
    'Snd10Days': (
        40,
        'day',
        ['snd'],
        dSndLE10Days_,
        ['year'],
        (1, None, (), {}),
        dict(
            name='days with 0-10cm snow-depth',
            units='days',
            attrU={'CLIMI': 'total days with 0 < snd <= 10 cm'}
            ),
        'p'),
    'Snd20Days': (
        41,
        'day',
        ['snd'],
        dSndGT10LE20Days_,
        ['year'],
        (1, None, (), {}),
        dict(
            name='days with 10-20cm snow-depth',
            units='days',
            attrU={'CLIMI': 'total days with 10 < snd <= 20 cm'}
            ),
       'p'),
    'SNWmax': (
        42,
        'day',
        ['snw'],
        'MAX',
        ['year'],
        (0, None),
        dict(name='maximum surface snow amount'),
        'p'),
    'TAS': (
        43,
        'mon',
        ['tas'],
        None,
        ['year'],
        (0, None),
        None,
        't'),
    'DTR': (
        44,
        'mon',
        ['tasmax', 'tasmin'],
        mDTR_,
        ['month'],
        (1, None),
        dict(
            name='daily temperature range',
            attrU={'CLIMI': 'tasmax - tasmin'}
            ),
        't'),
    'Rho925': (
        45,
        'day',
        ['ta925', 'hus925', 'ps'],
        None,
        ['month'],
        (4, ['ta925', 'hus925'], (92500.,)),
        dict(
            name='air density at 925 hPa',
            units='kg m-3',
            attrU={'CLIMI': 'with grid cells where ps < 975 hPa masked if ps '
                            'is available'}
            ),
       'w'),
    'RhoS': (
        46,
        'day',
        ['tas', 'huss', 'ps'],
        None,
        ['month'],
        (5, None),
        dict(name='near-surface air density', units='kg m-3'),
        'w'),
    'SuperCooledPR': (
        47,
        'day',
        ['tas', 'pr', 'ps',
         'ta925', 'ta850', 'ta700',
         'hus925', 'hus850', 'hus700'],
        dFreezRainDays_,
        ['year'],
        (2, None, 'y', {}, {}),
        dict(
            name='day with super cooled precipitation',
            units='days'
            ),
        'p'),
    'Snc25Days': (
        48,
        'day',
        ['snc'],
        dSncDays_,
        ['year'],
        (1, None, (), dict(thr=25)),
        dict(
            name='days with snow cover',
            units='days',
            attrU={'CLIMI': 'total days with snc > 25%'}
            ),
        'p'),
    'R5OScw': (
        49,
        'day',
        ['pr', 'tas', 'snc', 'snw'],
        dRainOnSnow_,
        ['year'],
        (6, ['dprrn', 'snc', 'snw'], (), {}),
        dict(
            name='days with rain on snow',
            units='days',
            attrU={'CLIMI': 'total days with snc > 25% and daily prrn > 5 mm '
                            'and snw > 3 mm'}
            ),
        'p'),
    'R1OScw': (
        50,
        'day',
        ['pr', 'tas', 'snc', 'snw'],
        dRainOnSnow_,
        ['year'],
        (6, ['dprrn', 'snc', 'snw'], (), dict(thr_r=1.)),
        dict(
            name='days with rain on snow',
            units='days',
            attrU={'CLIMI': 'total days with snc > 25% and daily prrn > 1 mm '
                            'and snw > 3 mm'}
            ),
        'p'),
    'R5OSc': (
        51,
        'day',
        ['pr', 'tas', 'snc'],
        dRainOnSnow_,
        ['year'],
        (1, ['dprrn', 'snc'], (), {}),
        dict(
            name='days with rain on snow',
            units='days',
            attrU={'CLIMI': 'total days with snc > 25% and daily prrn > 5 mm'}
            ),
        'p'),
    'R1OSc': (
        52,
        'day',
        ['pr', 'tas', 'snc'],
        dRainOnSnow_,
        ['year'],
        (1, ['dprrn', 'snc'], (), dict(thr_r=1.)),
        dict(
            name='days with rain on snow',
            units='days',
            attrU={'CLIMI': 'total days with snc > 25% and daily prrn > 1 mm'}
            ),
        'p'),
    'PRSNmax': (
        53,
        'day',
        ['prsn'],
        'MAX',
        ['year'],
        (0, None),
        dict(name='maximum snowfall flux'),
        'p'),
    'CalmDays': (
        54,
        'day',
        ['sfcWind', 'uas', 'vas'],
        dCalmDays_,
        ['season', 'year'],
        (1, ['sfcWind'], (), {}),
        dict(
            name='calm days',
            units='days',
            attrU={'CLIMI': 'total days with sfcWind < 2 m s-1'}
            ),
        'w'),
    'ConCalmDays': (
        55,
        'day',
        ['sfcWind', 'uas', 'vas'],
        dConCalmDays_,
        ['year'],
        (2, ['sfcWind'], 'y', {}, {}),
        dict(
            name='consecutive calm days',
            units='days',
            attrU={'CLIMI': 'longest calm spell defined as consecutive days '
                            'with sfcWind < 2 m s-1'}
            ),
        'wc'),
    'Wind975': (
        56,
        'day',
        ['ua975', 'va975', 'ps'],
        None,
        ['month', 'season', 'year'],
        (0, ['Wind975']),
        dict(
            name='wind speed at 975 hPa',
            attrU={'CLIMI': 'with grid cells where ps < 975 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'Wind975toSfc': (
        57,
        'day',
        ['ua975', 'va975', 'sfcWind', 'ps', 'uas', 'vas'],
        None,
        ['season', 'year'],
        (7, ['Wind975', 'sfcWind']),
        dict(
            name='ratio of wind speed at 975 hPa to surface wind speed',
            units=1,
            attrU={'CLIMI': 'with grid cells where ps < 975 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'ColdRainDays': (
        58,
        'day',
        ['pr', 'tas'],
        dColdRainDays_,
        ['year'],
        (1, None, (), {}),
        dict(
            name='days with cold rain',
            units='days',
            attrU={'CLIMI': 'total days with daily prrn > 0 mm and '
                            '0.58 degree C < tas < 2 degree C'}
            ),
        'p'),
    'ColdRainGT10Days': (
        59,
        'day',
        ['pr', 'tas'],
        dColdRainDays_,
        ['year'],
        (1, None, (), dict(thr_pr=10)),
        dict(
            name='days with cold rain',
            units='days',
            attrU={'CLIMI': 'total days with daily prrn > 10 mm and '
                            '0.58 degree C < tas < 2 degree C'}
            ),
        'p'),
    'ColdRainGT20Days': (
        60,
        'day',
        ['pr', 'tas'],
        dColdRainDays_,
        ['year'],
        (1, None, (), dict(thr_pr=20)),
        dict(
            name='days with cold rain',
            units='days',
            attrU={'CLIMI': 'total days with daily prrn > 20 mm and '
                            '0.58 degree C < tas < 2 degree C'}
            ),
        'p'),
    'WarmSnowDays': (
        61,
        'day',
        ['pr', 'tas'],
        dWarmSnowDays_,
        ['year'],
        (1, None, (), {}),
        dict(
            name='days with warm snow',
            units='days',
            attrU={'CLIMI': 'total days with daily prsn > 0 mm and '
                            '-2 degree C < tas < 0.58 degree C'}
            ),
        'p'),
    'WarmSnowGT10Days': (
        62,
        'day',
        ['pr', 'tas'],
        dWarmSnowDays_,
        ['year'],
        (1, None, (), dict(thr_pr=10)),
        dict(
            name='days with warm snow',
            units='days',
            attrU={'CLIMI': 'total days with daily prsn > 10 mm and '
                            '-2 degree C < tas < 0.58 degree C'}
            ),
        'p'),
    'WarmSnowGT20Days': (
        63,
        'day',
        ['pr', 'tas'],
        dWarmSnowDays_,
        ['year'],
        (1, None, (), dict(thr_pr=20)),
        dict(
            name='days with warm snow',
            units='days',
            attrU={'CLIMI': 'total days with daily prsn > 20 mm and '
                            '-2 degree C < tas < 0.58 degree C'}
            ),
        'p'),
    'ColdPRRNdays': (
        64,
        'day',
        ['pr', 'tas', 'prsn'],
        dColdPRRNdays_,
        ['year'],
        (1, None, (), {}),
        dict(
            name='days with cold rain',
            units='days',
            attrU={'CLIMI': 'total days with daily prrn > 0 mm and '
                            'tas < 2 degree C'}
            ),
        'p'),
    'ColdPRRNgt10Days': (
        65,
        'day',
        ['pr', 'tas', 'prsn'],
        dColdPRRNdays_,
        ['year'],
        (1, None, (), dict(thr_pr=10)),
        dict(
            name='days with cold rain',
            units='days',
            attrU={'CLIMI': 'total days with daily prrn > 10 mm and '
                            'tas < 2 degree C'}
            ),
        'p'),
    'ColdPRRNgt20Days': (
        66,
        'day',
        ['pr', 'tas', 'prsn'],
        dColdPRRNdays_,
        ['year'],
        (1, None, (), dict(thr_pr=20)),
        dict(
            name='days with cold rain',
            units='days',
            attrU={'CLIMI': 'total days with daily prrn > 20 mm and '
                            'tas < 2 degree C'}
            ),
        'p'),
    'WarmPRSNdays': (
        67,
        'day',
        ['tas', 'prsn'],
        dWarmPRSNdays_,
        ['year'],
        (1, None, (), {}),
        dict(
            name='days with warm snow',
            units='days',
            attrU={'CLIMI': 'total days with daily prsn > 0 mm and '
                            'tas > -2 degree C'}
            ),
        'p'),
    'WarmPRSNgt10Days': (
        68,
        'day',
        ['tas', 'prsn'],
        dWarmPRSNdays_,
        ['year'],
        (1, None, (), dict(thr_pr=10)),
        dict(
            name='days with warm snow',
            units='days',
            attrU={'CLIMI': 'total days with daily prsn > 10 mm and '
                            'tas > -2 degree C'}
            ),
        'p'),
    'WarmPRSNgt20Days': (
        69,
        'day',
        ['tas', 'prsn'],
        dWarmPRSNdays_,
        ['year'],
        (1, None, (), dict(thr_pr=20)),
        dict(
            name='days with warm snow',
            units='days',
            attrU={'CLIMI': 'total days with daily prsn > 20 mm and '
                            'tas > -2 degree C'}
            ),
        'p'),
    'SST': (
        70,
        'mon',
        ['tos'],
        None,
        ['season', 'year'],
        (0, None),
        None,
        't'),
    'SIC': (
        71,
        'mon',
        ['sic'],
        None,
        ['season', 'year'],
        (0, None),
        None,
        'tp'),
    'Rho975': (
        72,
        'day',
        ['ta975', 'hus975', 'ps'],
        None,
        ['month'],
        (4, ['ta975', 'hus975'], (97500.,)),
        dict(
            name='air density at 975 hPa',
            units='kg m-3',
            attrU={'CLIMI': 'with grid cells where ps < 975 hPa masked if ps '
                            'is available'}),
        'w'),
    'h6SuperCooledPR': (
        73,
        '6hr',
        ['tas', 'pr', 'ps',
         'ta925', 'ta850', 'ta700',
         'hus925', 'hus850', 'hus700'],
        dFreezRainDays_,
        ['year'],
        (2, None, 'y', {}, {}),
        dict(
            name='number of super cooled precipitation events',
            units=1,
            attrU={'CLIMI': 'interval of 6hr'}
            ),
        'p'),
    'Wind925': (
        74,
        'day',
        ['ua925', 'va925', 'ps'],
        None,
        ['month', 'season', 'year'],
        (0, ['Wind925']),
        dict(
            name='wind speed at 925 hPa',
            attrU={'CLIMI': 'with grid cells where ps < 925 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'Wind925toSfc': (
        75,
        'day',
        ['ua925', 'va925', 'sfcWind', 'ps', 'uas', 'vas'],
        None,
        ['season', 'year'],
        (7, ['Wind925', 'sfcWind']),
        dict(
            name='ratio of wind speed at 925 hPa to surface wind speed',
            units=1,
            attrU={'CLIMI': 'with grid cells where ps < 925 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'FirstDayWithoutFrost': (
        76,
        'day',
        ['tasmin'],
        dFirstDayWithoutFrost_,
        ['year'],
        (1, None, 'yd', {}, {}),
        dict(
            name='first day without frost',
            units=1,
            attrU={'CLIMI': 'first day with tasmin > 0 degree C'}
            ),
        't'),
    'CalmDays925': (
        77,
        'day',
        ['ua925', 'va925', 'ps'],
        dCalmDays_,
        ['year'],
        (1, ['Wind925'], (), {}),
        dict(
            name='calm days at 925 hPa',
            units='days',
            attrU={'CLIMI': 'total days with wind speed at 925 hPa '
                            '< 2 m s-1; '
                            'with grid cells where ps < 925 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'ConCalmDays925': (
        78,
        'day',
        ['ua925', 'va925', 'ps'],
        dConCalmDays_,
        ['year'],
        (2, ['Wind925'], 'y', {}, {}),
        dict(
            name='consecutive calm days at 925 hPa',
            units='days',
            attrU={'CLIMI': 'longest calm spell defined as consecutive days '
                            'with wind speed at 925 hPa < 2 m s-1; '
                            'with grid cells where ps < 925 hPa masked if ps '
                            'is available'}
            ),
        'wc'),
    'CalmDays975': (
        79,
        'day',
        ['ua975', 'va975', 'ps'],
        dCalmDays_,
        ['year'],
        (1, ['Wind975'], (), {}),
        dict(
            name='calm days at 975 hPa',
            units='days',
            attrU={'CLIMI': 'total days with wind speed at 975 hPa '
                            '< 2 m s-1; '
                            'with grid cells where ps < 975 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'ConCalmDays975': (
        80,
        'day',
        ['ua975', 'va975', 'ps'],
        dConCalmDays_,
        ['year'],
        (2, ['Wind975'], 'y', {}, {}),
        dict(
            name='consecutive calm days at 975 hPa',
            units='days',
            attrU={'CLIMI': 'longest calm spell defined as consecutive days '
                            'with wind speed at 975 hPa < 2 m s-1; '
                            'with grid cells where ps < 975 hPa masked if ps '
                            'is available'}
            ),
        'wc'),
    'MinusDays': (
        81,
        'day',
        ['tas'],
        dMinusDays_,
        ['year'],
        (1, None),
        dict(
            name='minus days',
            units='days',
            attrU={'CLIMI': 'total days with tas < 0 degree C'}
            ),
        't'),
    'FreezingDays': (
        82,
        'day',
        ['tasmax'],
        dFreezingDays_,
        ['year'],
        (1, None),
        dict(
            name='freezing days',
            units='days',
            attrU={'CLIMI': 'total days with tasmax < 0 degree C'}
            ),
       't'),
    'ColdRainWarmSnowDays': (
        83,
        'day',
        ['pr', 'tas'],
        dColdRainWarmSnowDays_,
        ['year'],
        (1, None),
        dict(
            name='days with cold rain or warm snow',
            units='days',
            attrU={'CLIMI': 'total days with daily pr > 0 mm and '
                            '-2 degree C < tas < 2 degree C'}
            ),
        'p'),
    'Wind50m': (
        84,
        'day',
        ['ua50m', 'va50m'],
        None,
        ['month', 'season', 'year'],
        (0, ['Wind50m']),
        dict(name='wind speed at 50 m height'),
        'w'),
    'Rho50m': (
        85,
        'day',
        ['ta50m', 'hus50m', 'p50m'],
        None,
        ['month'],
        (5, None),
        dict(name='air density at 50 m height', units='kg m-3'),
        'w'),
    'Wind100m': (
        86,
        'day',
        ['ua100m', 'va100m'],
        None,
        ['month', 'season', 'year'],
        (0, ['Wind100m']),
        dict(name='wind speed at 100 m height'),
        'w'),
    'Rho100m': (
        87,
        'day',
        ['ta100m', 'hus100m', 'p100m'],
        None,
        ['month'],
        (5, None),
        dict(name='air density at 100 m height', units='kg m-3'),
        'w'),
    'Wind200m': (
        88,
        'day',
        ['ua200m', 'va200m'],
        None,
        ['month', 'season', 'year'],
        (0, ['Wind200m']),
        dict(name='wind speed at 200 m height'),
        'w'),
    'Rho200m': (
        89,
        'day',
        ['ta200m', 'hus200m', 'p200m'],
        None,
        ['month'],
        (5, None),
        dict(name='air density at 200 m height', units='kg m-3'),
        'w'),
    'Wind950': (
        90,
        'day',
        ['ua950', 'va950', 'ps'],
        None,
        ['month', 'season', 'year'],
        (0, ['Wind950']),
        dict(
            name='wind speed at 950 hPa',
            attrU={'CLIMI': 'with grid cells where ps < 950 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'Rho950': (
        91,
        'day',
        ['ta950', 'hus950', 'ps'],
        None,
        ['month'],
        (4, ['ta950', 'hus950'], (95000.,)),
        dict(
            name='air density at 950 hPa',
            units='kg m-3',
            attrU={'CLIMI': 'with grid cells where ps < 950 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'Wind950toSfc': (
        92,
        'day',
        ['ua950', 'va950', 'sfcWind', 'ps', 'uas', 'vas'],
        None,
        ['season', 'year'],
        (7, ['Wind950', 'sfcWind']),
        dict(
            name='ratio of wind speed at 950 hPa to surface wind speed',
            units=1,
            attrU={'CLIMI': 'with grid cells where ps < 950 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'Wind900': (
        93,
        'day',
        ['ua900', 'va900', 'ps'],
        None,
        ['month', 'season', 'year'],
        (0, ['Wind900']),
        dict(
            name='wind speed at 900 hPa',
            attrU={'CLIMI': 'with grid cells where ps < 900 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'Rho900': (
        94,
        'day',
        ['ta900', 'hus900', 'ps'],
        None,
        ['month'],
        (4, ['ta900', 'hus900'], (90000.,)),
        dict(
            name='air density at 900 hPa', units='kg m-3',
            attrU={'CLIMI': 'with grid cells where ps < 900 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'Wind900toSfc': (
        95,
        'day',
        ['ua900', 'va900', 'sfcWind', 'ps', 'uas', 'vas'],
        None, ['season', 'year'],
        (7, ['Wind900', 'sfcWind']),
        dict(
            name='ratio of wind speed at 900 hPa to surface wind speed',
            units=1,
            attrU={'CLIMI': 'with grid cells where ps < 900 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'CalmDays50m': (
        96,
        'day',
        ['ua50m', 'va50m'],
        dCalmDays_,
        ['year'],
        (1, ['Wind50m'], (), {}),
        dict(
            name='calm days at 50 m height',
            units='days',
            attrU={'CLIMI': 'total days with wind speed at 50 m height '
                            '< 2 m s-1'}
            ),
        'w'),
    'ConCalmDays50m': (
        97,
        'day',
        ['ua50m', 'va50m'],
        dConCalmDays_,
        ['year'],
        (2, ['Wind50m'], 'y', {}, {}),
        dict(
            name='consecutive calm days at 50 m height',
            units='days',
            attrU={'CLIMI': 'longest calm spell defined as consecutive days '
                            'with wind speed at 50 m height < 2 m s-1'}
            ),
        'wc'),
    'CalmDays100m': (
        98,
        'day',
        ['ua100m', 'va100m'],
        dCalmDays_,
        ['year'],
        (1, ['Wind100m'], (), {}),
        dict(
            name='calm days at 100 m height',
            units='days',
            attrU={'CLIMI': 'total days with wind speed at 100 m height '
                            '< 2 m s-1'}
            ),
        'w'),
    'ConCalmDays100m': (
        99,
        'day',
        ['ua100m', 'va100m'],
        dConCalmDays_,
        ['year'],
        (2, ['Wind100m'], 'y', {}, {}),
        dict(
            name='consecutive calm days at 100 m height',
            units='days',
            attrU={'CLIMI': 'longest calm spell defined as consecutive days '
                            'with wind speed at 100 m height < 2 m s-1'}
            ),
        'wc'),
    'CalmDays200m': (
        100,
        'day',
        ['ua200m', 'va200m'],
        dCalmDays_,
        ['year'],
        (1, ['Wind200m'], (), {}),
        dict(
            name='calm days at 200 m height',
            units='days',
            attrU={'CLIMI': 'total days with wind speed at 200 m height '
                            '< 2 m s-1'}
            ),
        'w'),
    'ConCalmDays200m': (
        101,
        'day',
        ['ua200m', 'va200m'],
        dConCalmDays_,
        ['year'],
        (2, ['Wind200m'], 'y', {}, {}),
        dict(
            name='consecutive calm days at 200 m height',
            units='days',
            attrU={'CLIMI': 'longest calm spell defined as consecutive days '
                            'with wind speed at 200 m height < 2 m s-1'}
            ),
        'wc'),
    'CalmDays950': (
        102,
        'day',
        ['ua950', 'va950', 'ps'],
        dCalmDays_,
        ['year'],
        (1, ['Wind950'], (), {}),
        dict(
            name='calm days at 950 hPa',
            units='days',
            attrU={'CLIMI': 'total days with wind speed at 950 hPa '
                            '< 2 m s-1; '
                            'with grid cells where ps < 950 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'ConCalmDays950': (
        103,
        'day',
        ['ua950', 'va950', 'ps'],
        dConCalmDays_,
        ['year'],
        (2, ['Wind950'], 'y', {}, {}),
        dict(
            name='consecutive calm days at 950 hPa',
            units='days',
            attrU={'CLIMI': 'longest calm spell defined as consecutive days '
                            'with wind speed at 950 hPa < 2 m s-1; '
                            'with grid cells where ps < 950 hPa masked if ps '
                            'is available'}
            ),
        'wc'),
    'CalmDays900': (
        104,
        'day',
        ['ua900', 'va900', 'ps'],
        dCalmDays_,
        ['year'],
        (1, ['Wind900'], (), {}),
        dict(
            name='calm days at 900 hPa',
            units='days',
            attrU={'CLIMI': 'total days with wind speed at 900 hPa '
                            '< 2 m s-1; '
                            'with grid cells where ps < 900 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'ConCalmDays900': (
        105,
        'day',
        ['ua900', 'va900', 'ps'],
        dConCalmDays_, ['year'],
        (2, ['Wind900'], 'y', {}, {}),
        dict(
            name='consecutive calm days at 900 hPa',
            units='days',
            attrU={'CLIMI': 'longest calm spell defined as consecutive days '
                            'with wind speed at 900 hPa < 2 m s-1; '
                            'with grid cells where ps < 900 hPa masked if ps '
                            'is available'}
            ),
        'wc'),
    'h3SfcWind': (
        106,
        '3hr',
        ['sfcWind', 'uas', 'vas'],
        None,
        ['hour-month', 'hour-season', 'hour'],
        (0, ['sfcWind']),
        None,
        'w'),
    'h3RhoS': (
        107,
        '3hr',
        ['tas', 'huss', 'ps'],
        None,
        ['hour-month'],
        (5, None),
        dict(name='near-surface air density', units='kg m-3'),
        'w'),
    'h3Wind975': (
        108,
        '3hr',
        ['ua975', 'va975', 'ps'],
        None,
        ['hour-month', 'hour-season', 'hour'],
        (8, ['ua975','va975'], (97500.,)),
        dict(
            name='wind speed at 975 hPa',
            attrU={'CLIMI': 'with grid cells where ps < 975 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'h3Rho975': (
        109,
        '3hr',
        ['ta975', 'hus975', 'ps'],
        None,
        ['hour-month'],
        (4, ['ta975', 'hus975'], (97500.,)),
        dict(
            name='air density at 975 hPa',
            units='kg m-3',
            attrU={'CLIMI': 'with grid cells where ps < 975 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'h3Wind950': (
        110,
        '3hr',
        ['ua950', 'va950', 'ps'],
        None,
        ['hour-month', 'hour-season', 'hour'],
        (8, ['ua950','va950'], (95000.,)),
        dict(
            name='wind speed at 950 hPa',
            attrU={'CLIMI': 'with grid cells where ps < 950 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'h3Rho950': (
        111,
        '3hr',
        ['ta950', 'hus950', 'ps'],
        None,
        ['hour-month'],
        (4, ['ta950', 'hus950'], (95000.,)),
        dict(
            name='air density at 950 hPa',
            units='kg m-3',
            attrU={'CLIMI': 'with grid cells where ps < 950 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'h3Wind925': (
        112,
        '3hr',
        ['ua925', 'va925', 'ps'],
        None,
        ['hour-month', 'hour-season', 'hour'],
        (8, ['ua925','va925'], (92500.,)),
        dict(
            name='wind speed at 925 hPa',
            attrU={'CLIMI': 'with grid cells where ps < 925 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'h3Rho925': (
        113,
        '3hr',
        ['ta925', 'hus925', 'ps'],
        None,
        ['hour-month'],
        (4, ['ta925', 'hus925'], (92500.,)),
        dict(
            name='air density at 925 hPa',
            units='kg m-3',
            attrU={'CLIMI': 'with grid cells where ps < 975 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'h3Wind900': (
        114,
        '3hr',
        ['ua900', 'va900', 'ps'],
        None,
        ['hour-month', 'hour-season', 'hour'],
        (8, ['ua900','va900'], (90000.,)),
        dict(
            name='wind speed at 900 hPa',
            attrU={'CLIMI': 'with grid cells where ps < 900 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'h3Rho900': (
        115,
        '3hr',
        ['ta900', 'hus900', 'ps'],
        None,
        ['hour-month'],
        (4, ['ta900', 'hus900'], (90000.,)),
        dict(
            name='air density at 900 hPa',
            units='kg m-3',
            attrU={'CLIMI': 'with grid cells where ps < 900 hPa masked if ps '
                            'is available'}
            ),
        'w'),
    'h3Wind50m': (
        116,
        '3hr',
        ['ua50m', 'va50m'],
        None,
        ['hour-month', 'hour-season', 'hour'],
        (0, ['Wind50m']),
        dict(
            name='wind speed at 50 m height'),
        'w'),
    'h3Rho50m': (
        117,
        '3hr',
        ['ta50m', 'hus50m', 'p50m'],
        None,
        ['hour-month'],
        (5, None),
        dict(name='air density at 50 m height', units='kg m-3'),
        'w'),
    'h3Wind100m': (
        118,
        '3hr',
        ['ua100m', 'va100m'],
        None,
        ['hour-month', 'hour-season', 'hour'],
        (0, ['Wind100m']),
        dict(name='wind speed at 100 m height'),
        'w'),
    'h3Rho100m': (
        119,
        '3hr',
        ['ta100m', 'hus100m', 'p100m'],
        None,
        ['hour-month'],
        (5, None),
        dict(name='air density at 100 m height', units='kg m-3'),
        'w'),
    'h3Wind200m': (
        120,
        '3hr',
        ['ua200m', 'va200m'],
        None,
        ['hour-month', 'hour-season', 'hour'],
        (0, ['Wind200m']),
        dict(name='wind speed at 200 m height'),
        'w'),
    'h3Rho200m': (
        121,
        '3hr',
        ['ta200m', 'hus200m', 'p200m'],
        None,
        ['hour-month'],
        (5, None),
        dict(name='air density at 200 m height', units='kg m-3'),
        'w')
    }


v__ = {
    'sfcWind': (
        ws_cube,
        ['uas', 'vas'],
        ((), {}),
        None),
    'Wind50m': (
        ws_cube,
        ['ua50m', 'va50m'],
        ((), {}),
        None),
    'Wind100m': (
        ws_cube,
        ['ua100m', 'va100m'],
        ((), {}),
        None),
    'Wind200m': (
        ws_cube,
        ['ua200m', 'va200m'],
        ((), {}),
        None),
    'Wind900': (
        ws_cube,
        ['ua900', 'va900'],
        ((), {}),
        mask_cube,
        90000.),
    'Wind925': (
        ws_cube,
        ['ua925', 'va925'],
        ((), {}),
        mask_cube,
        92500.),
    'Wind950': (
        ws_cube,
        ['ua950', 'va950'],
        ((), {}),
        mask_cube,
        95000.),
    'Wind975': (
        ws_cube,
        ['ua975', 'va975'],
        ((), {}),
        mask_cube,
        97500.),
    'dprrn': (
        dPRRN_fr_PR_T_,
        ['pr', 'tas'],
        ((), {}),
        None),
    'mrro_amjjas': (
        extract_season_cube,
        ['mrro'],
        (('amjjas',), {}),
        None),
    'pMINUSe': (
        'c_pr.copy(c_pr.data - c_evspsbl.data)',),
    'prrn': (
        'c_pr.copy(c_pr.data - c_prsn.data)',)
    }
