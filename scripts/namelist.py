from datetime import datetime
import numpy as np
import os


# The variable "namelist" defines the information content on which the
# analyses are based. 
#
# "namelist" is a two-item list:
    # 1. The first item is a list of starting dates (datetime format)
    # 2. The second item is a list of contributors
        # Each contributor may or may have not contributed forecasts for 
        # each season. A "contributor" item is a 1 + n item list where n is
        # the number of starting dates. The first item in contributor is its 
        # consistent name


# Information on the dates. Each raw corresponds to one call / outlook / exercise
#                     Start date             End of forecast period   Start - end verification period                    Start - end forecast data availability period
namelistOutlooks = [ [datetime(2017, 12, 1), datetime(2018, 2, 28)  , datetime(2018, 2, 1)     , datetime(2018, 2, 28),  datetime(2018, 2, 1), datetime(2018, 2, 28) ],
                     [datetime(2018, 12, 1), datetime(2019, 2, 28)  , datetime(2019, 2, 1)     , datetime(2019, 2, 28),  datetime(2018, 12,1), datetime(2019, 2, 28) ],
                     [datetime(2019, 12, 1), datetime(2020, 2, 28)  , datetime(2020, 2, 1)     , datetime(2020, 2, 28),  datetime(2019, 12,1), datetime(2020, 2, 28) ],
                     [datetime(2020, 12, 1), datetime(2021, 2, 28)  , datetime(2021, 2, 1)     , datetime(2021, 2, 28),  datetime(2020, 12,1), datetime(2021, 2, 28) ],
                     [datetime(2021, 12, 1), datetime(2022, 2, 28)  , datetime(2022, 2, 1)     , datetime(2022, 2, 28),  datetime(2021, 12,1), datetime(2022, 2, 28) ],
                     [datetime(2022, 12, 1), datetime(2023, 2, 28)  , datetime(2023, 2, 1)     , datetime(2023, 2, 28),  datetime(2022, 12,1), datetime(2023, 2, 28) ],
                     [datetime(2023, 12, 1), datetime(2024, 2, 28)  , datetime(2024, 2, 1)     , datetime(2024, 2, 28),  datetime(2023, 12,1), datetime(2024, 2, 28) ]
                   ]

#                                       2017 - 2018               2018-2019             2019-2020                 2020-2021              2021-2022             2022-2023         2023-2024
#         Short name                nb. forecasts for 
#                                   SIA, rSIA, SIC, SIV      SIA, rSIA, SIC, SIV    SIA, rSIA, SIC, SIV   SIA, rSIA, SIC, SIV    SIA, rSIA, SIC, SIV    SIA, rSIA, SIC, SIV     SIA, rSIA, SIC, SIV


namelistContributions = [
              ["AWI-SDAP",      	[0,   0,   0,    0  ],  [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   1,    0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ],  [0,   0,    0,   0  ],      "s"  ], \
              ["AWI-CPS",           [0,   0,   0,    0  ],  [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ],  [30,   30,   30,  30  ],    "s"  ], \
              ["barreira",         	[0,   0,   0,    0  ],  [1,   1,   1,    0  ], [1,   1,   1,    0  ], [1,   1,   1,    0  ], [3,   3,    3,   0  ], [2,   2,    2,   0  ],  [0,   0,    0,   0  ],      "s"  ], \
              ["BSC"        ,      	[0,   0,   0,    0  ],  [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [10,  10,  10,   0  ], [10,  10,   10,  0  ],  [50,  50,   50,  0  ],      "d"  ], \
              ["CanSIPSv2",        	[0,   0,    0,   0  ],  [0,   0,    0,   0  ], [0,   0,    0,   0  ], [20,  20,   0,   0  ], [20,  20,   0,   0  ], [20,  20,   0,   0  ],  [20,  20,   0,   0  ],      "d"  ], \
              ["climatology"  ,    	[30,  0 ,  30,   0],    [30,  0 ,  30,   0  ], [30,  0 ,  30,   0  ], [30,  0 ,  30,   0  ], [30,  0 ,  30,   0  ], [30,  0 ,  30,   0  ],  [30,  0 ,  30,   0  ],      "b"  ], \
              ["cmcc",            	[0,   0,   0,    0  ],  [50,  0,   0,    0  ], [0,   0,   0,    0  ], [50,  50,  50,   0  ], [50,  50,   50,  0  ], [50,  50,   50,  0  ],  [50,  50,   50,  0  ],      "d"  ], \
              ["CNRM",     	        [0,   0,   0,    0  ],  [0,   0,   0,    0  ], [51,  51,  51,   0  ], [51,  51,  51,   0  ], [51,  51,   51,  51 ], [0,   0,   0,    0  ],  [0,   0,   0,    0  ],      "d"  ], \
              ["ecmwf",           	[50,  50,   0,   0],    [50,  50,   0,   0  ], [51,  51,   0,   0  ], [51,  51,   0,   0  ], [51,  51,    0,   0 ], [51,  51,    0,   0 ],  [0,  0,    0,   0 ],        "d"  ], \
              ["emc",             	[15,  15,  0,    0  ],  [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,    0,   0  ], [10,  10,  10,    0 ],  [10,  10,  10,    0  ],     "d"  ], \
              ["FIO-ESM",         	[1,   0,    0,   0  ], 	[1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ], [10,   0,    0,   0 ],  [0,   0,    0,   0 ],       "d"  ], \
              ["Gateway",         	[1,   0,    0,   0  ], 	[0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ],  [0,   0,    0,   0  ],      "s"  ], \
              ["gfdl"    ,        	[0,   0,   0,    0  ],  [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [30,  30,   30,  0  ], [30,  30,   30,  0  ],  [30,  30,   30,  0  ],      "d"  ], \
              ["IOCAS-SIPNet",      [0,   0,   0,    0  ],  [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ],  [1,   0,   0,    0  ],      "s"  ], \
              ["Lamont",          	[1,   1,    0,   0  ],  [1,   1,    1,   0  ], [1,   1,    1,   0  ], [1,   1,    1,   0  ], [1,   1,    1,   0  ], [1,   1,    1,   0  ],  [1,   1,    1,   0  ],      "s"  ], \
              ["Meier-NSIDC",      	[0,   0,   0,    0  ],  [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [1,   0,   0,    0  ], [1,   0,   0,    0  ],  [1,   0,   0,    0  ],      "s"  ], \
              ["MetOffice",        	[42,  42,  42,   0],    [42,  42,  42,   0  ], [42,  42,  42,   0  ], [42,  42,  42,   0  ], [42,  42,  42,   0  ], [42,  42,  42,   0  ],  [42,  42,  42,   0  ],      "d"  ], \
              ["Modified-CanSIPS",  [20,  0,    0,   0  ],  [20,  20,   0,   0  ], [20,  20,   0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ],  [0,   0,    0,   0  ],      "d"  ], \
              ["mpas-cesm",        	[2,   2,    0,   0  ], 	[0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ],  [0,   0,    0,   0  ],      "d"  ], \
              ["nasa-gmao",        	[10,  10,  10,   0 ],   [10,  0,   10,   0, ], [10,  0,   10,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ],  [0,   0,    0,   0  ],      "d"  ], \
              ["NASA-GSFC",       	[1,   0,    0,   0  ],  [1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ],  [0,   0,    0,   0  ],      "s"  ], \
              ["NicoSun",          	[1,   0,    3,   0  ],  [3,   0,    3,   0  ], [3,   0,    3,   0  ], [3,   0,    3,   3  ], [3,   0,    3,   3  ], [3,   0,    3,   3  ],  [3,   0,    3,   3  ],      "s"  ], \
              ["nrl",              	[0,   0,    0,   0  ],  [0,   0,   10,  0   ],  [0,   0,    0,   0 ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ],  [0,   0,    0,   0  ],      "d"  ], \
              ["SINTEX-F2",        	[0,   0,   0,    0  ],  [0,   0,   0,    0  ], [0,   0,   0,    0  ], [12,  12,  0,   0   ], [24,  24,   0,   0  ], [24,  24,   0,   0  ],  [0,   0,    0,   0  ],      "d"  ], \
              ["SYSU-SML-KNN" ,     [0,   0,   0,    0  ], 	[0,   0,   0,    0  ], [0,   0,   0,    0  ], [1,   1,   1,   0   ], [1,   1,    1,   0  ], [1,   0,    1,   0  ],  [1,   0,    0,   0  ],      "s"  ], \
              ["SYSU-SML-ConvLSTM", [0,   0,   0, 0 ], 	    [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,   0   ], [0,   0,    0,   0  ], [1,   1,    0,   0  ],  [1,   0,    1,   0  ],      "s"  ], \
              ["SYSU-SML-MLM",      [0,   0,   0,    0  ],  [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,   0   ], [0,   0,    0,   0  ], [1,   0,    1,   0  ],  [1,   0,    0,   0  ],      "s"  ], \
              ["ucl",             	[10,  10,  10,   0  ],  [10,  10,  10,   0  ], [10,  10,  10,   0  ], [10,  10,  10,   10 ], [10,  10,  10,   10 ], [10,  10,  10,   10 ],  [10,  10,  10,   10 ],      "d"  ], \
              ["UW",               	[0,   0,   0,    0  ], 	[ 0,   0,   0,   0  ], [ 0,   0,   0,   0  ], [ 0,   0,   0,    0 ], [ 0,   0,   0,    0 ], [10,  0,    0,    0 ],  [0,  0,    0,    0 ],       "d"  ], \
           ]
    

#           Name                       Lon western  Lon eastern

# Don't change names! (used in scripts)
sectors = [
           ["Weddell",                 -60.0,       20.0  ],
           ["Indian",                  20.0 ,       90.0  ], 
           ["West Pacific",            90.0 ,       160.0 ],
           ["Ross",                    160.0,       -130.0],
           ["Amundsen-Bellingshausen", -130.0,      -60.0 ],
           ["Southern Ocean"         , -180.0,      180.0 ],
          ]
           
# def check_files():
    
#     diags      = ["total-area", "regional-area", "concentration", "volume"]
#     nDiags     = len(diags)
#     extensions = ["txt"       , "txt",           "nc",            "nc"    ]
#     folder     = ["txt",        "txt",           "netcdf",        "netcdf"]
#     for sub in namelist[1]:

    
#         print("Checking data availability for " + sub[0])
        
#         for ji, iniDate in enumerate(namelist[0]):
#             for jd in range(nDiags): # diag
            
                
#                 nFiles = sub[ji + 1][jd]
                     
#                 if nFiles > 0:
#                     for jFile in np.arange(1, nFiles + 1):
#                         folder = str(iniDate.year) + "-" + str(iniDate.year + 1) + "/" + extensions[jd] + "/"
                        
#                         file = "../data/" + "/" + folder + "/" + sub[0] + "_" + str(jFile).zfill(3) + "_" + diags[jd] + "." + extensions[jd]
#                         os.path.isfile(file)
#                         if os.path.isfile(file):
#                             pass
#                         else:
#                             print("   " + file + " not found")
                
            
            
        
        
# check_files()
