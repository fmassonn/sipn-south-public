#          Group               Forecast ID      Color for plots [string]     Comment
#          name [string]       list [list]
np.random.seed(10)
info = [ 
    [ "climatology",      range(1, 30  + 1),   [0.1, 0.1, 0.1],  "",        ],   \ 
    [ "AWI-CPS",          range(1, 30),        [0.0, 0.0, 1.0],  "$^d$",    ],   \
    [ "BSC",              range(1, 50 + 1),    [0.6, 0.8, 0.9],  "$^d$",    ],   \
    [ "cmcc",             range(1, 50 + 1),    "#6495ED",        "$^d$",    ],   \
#    [ "emc",              range(1, 10 + 1),    "#00008B",        "$^d$"     ],   \
    [ "gfdl",             range(1, 30 + 1),    "#5D3FD3",        "$^d$"     ],   \
    [ "IOCAS-SIPNet",     range(1, 5 + 1),     [1.0, 0.75, 0.0], "$^s$"     ], \
    [ "Lamont",           range(1, 3)    ,     [1.0, 0.78, 0.6], "$^{s,i}$" ], \
    [ "MetOffice",        range(1, 42 + 1),    [0.0, 1.0, 1.0],  "$^d$"     ], \
#    [ "NicoSun",          range(1, 3 + 1),     "#CD7F32",        "$^s$"     ], \
    [ "SYSU-SML-KNN",     [1],                 "#CC5500",        "$^s$"    ], \
    [ "SYSU-SML-MLM",     [1],                 "#E3963E",        "$^s$"    ], \
    [ "ucl",              range(1, 10 + 1),    "#A7C7E7",        "$^d$"    ], \
       ]


