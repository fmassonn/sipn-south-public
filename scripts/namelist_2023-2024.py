#          Group               Forecast ID    Tag       Area of ocean south of 60 S  
#          name [string]       list [list]
np.random.seed(10)
info = [ 
    [ "climatology",      range(1, 30  + 1),  ""  ],
     [ "AWI-CPS",          range(1, 30 + 1) ,                 "$^{d}$",  18.8333     ], \
#    [ "barreiraOrig",     range(1, 2  + 1),  "$^{s}$"  ], \
#    [ "barreiraNew",      range(1, 2  + 1),  "$^{s}$"  ], \
    [ "BSC",              range(1, 50 + 1),  "$^d$"    ], \
     [ "CanSIPSv2",        range(1, 20 + 1),  "$^d$",   20.24169    ], \
    [ "cmcc",             range(1, 50 + 1),  "$^d$"    ], \
#    [ "CNRM",             range(1, 51 + 1),  "$^d$"    ],  \
#    [ "ecmwf",            range(1, 51 + 1),  "$^d$"    ], \
    [ "emc",              range(1, 10 + 1),  "$^d$"    ], \
#    [ "FIO-ESM",          range(1, 10 + 1),  "$^d$"    ], \
    [ "gfdl",             range(1, 30 + 1),  "$^d$"    ], \
    ["IOCAS-SIPNet",      [1],               "$^s$"    ], \
   [ "Lamont",           [1],               "$^{s,i}$"], \
    [ "Meier-NSIDC",      [1],               "$^s$"    ], \
    [ "MetOffice",        range(1, 42 + 1),  "$^d$"    ], \
#    [ "NASA-GSFC",        [1],                "$^{s}$"  ],  \
    [ "NicoSun",          range(1, 3 + 1),   "$^s$"    ], \
#    [ "SINTEX-F2",       range(1, 24 + 1),   "$^d$"    ],  \
    [ "SYSU-SML-MLM",     [1],               "$^s$"    ], \
    [ "SYSU-SML-KNN",     [1],               "$^s$"    ], \
    [ "ucl",              range(1, 10 + 1),  "$^d$"    ], \
#    [ "UW" ,              range(1, 15 + 1),  "$^d$"    ], \
       ]


