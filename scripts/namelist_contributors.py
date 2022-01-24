import numpy as npimport os#                                    2017 - 2018             2018-2019         2019-2020                   2020-2021              2021-2022#         Short name           nb. forecasts for #                              SIA, rSIA, SIC, SIV    SIA, rSIA, SIC, SIV    SIA, rSIA, SIC, SIV   SIA, rSIA, SIC, SIV    SIA, rSIA, SIC, SIV# Years defining the dataset, corresponding to December of the relevant seasonyearStart = 2017yearEnd   = 2021info = [ ["nrl",              [6,   6,    6,   0  ], [10,  10,   10,  0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ]   ], \         ["NicoSun",          [1,   0,    3,   0  ], [3,   0,    3,   0  ], [3,   0,    3,   0  ], [3,   0,    3,   3  ], [3,   0,    3,   3  ]   ], \         ["FIO-ESM",          [1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ]   ], \         ["ecmwf",            [50,  50,   0,   0  ], [50,  50,   0,   0  ], [51,  51,   0,   0  ], [51,  51,   0,   0  ], [0,   0,    0,   0  ]   ], \         ["Gateway",          [1,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ]   ], \         ["mpas-cesm",        [2,   2,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ]   ], \         ["Modified-CanSIPS", [20,  0,    0,   0  ], [20,  20,   0,   0  ], [20,  20,   0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ]   ], \         ["CanSIPSv2",        [0,   0,    0,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ], [20,  20,   0,   0  ], [20,  20,   0,   0  ]   ], \         ["NASA-GSFC",        [1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ], [1,   0,    0,   0  ]   ], \         ["Lamont",           [1,   1,    1,   0  ], [1,   1,    1,   0  ], [1,   1,    1,   0  ], [1,   1,    1,   0  ], [1,   1,    1,   0  ]   ], \         ["nasa-gmao",        [10,  10,  10,   0  ], [10,  0,   10,   0, ], [10,  0,   10,   0  ], [0,   0,    0,   0  ], [0,   0,    0,   0  ]   ], \         ["MetOffice",        [42,  42,  42,   0  ], [42,  42,  42,   0  ], [42,  42,  42,   0  ], [42,  42,  42,   0  ], [42,  42,  42,   0  ]   ], \         ["ucl",              [10,  10,  10,   0  ], [10,  10,  10,   0  ], [10,  10,  10,   0  ], [10,  10,  10,   10 ], [10,  10,  10,   10 ]   ], \         ["emc",              [15,  15,  0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,    0,   0  ]   ], \         ["barreira",         [0,   0,   0,    0  ], [1,   1,   1,    0  ], [1,   1,   1,    0  ], [1,   1,   1,    0  ], [3,   3,    3,   0  ]   ], \         ["cmcc",             [0,   0,   0,    0  ], [50,  0,   0,    0  ], [0,   0,   0,    0  ], [50,  50,  50,   0  ], [50,  50,   50,  0  ]   ], \         ["CNRM",             [0,   0,   0,    0  ], [0,   0,   0,    0  ], [51,  51,  51,   0  ], [51,  51,  51,  51  ], [51,  51,   51,  51 ]   ], \         ["SINTEX-F2",        [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [12,  12,  0,   0   ], [24,  24,   0,   0  ]   ], \         ["SYSU"    ,         [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [1,   1,   0,   0   ], [1,   1,    0,   0  ]   ], \         ["gfdl"    ,         [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [30,  30,   30,  0  ]   ], \         ["Meier-NSIDC",      [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [0,   0,   0,    0  ], [1,   0,   0,    0  ]   ], \        ]            def check_files():        diags      = ["total-area", "regional-area", "concentration", "volume"]    nDiags     = len(diags)    extensions = ["txt"       , "txt",           "nc",            "nc"    ]    folder     = ["txt",        "txt",           "netcdf",        "netcdf"]    for i in info:        print("Checking data availability for " + i[0])                for j, y in enumerate(np.arange(yearStart, yearEnd + 1)):            for d in range(nDiags): # diag                            nFiles = i[j + 1][d]                                       if nFiles > 0:                    for jFile in np.arange(1, nFiles + 1):                        file = "../data/" + str(y) + "-" + str(y + 1) + "/" + folder[d] + "/" + i[0] + "_" + str(jFile).zfill(3) + "_" + diags[d] + "." + extensions[d]                        os.path.isfile(file)                        if os.path.isfile(file):                            pass                        else:                            print("   " + file + " not found")                                                        check_files()