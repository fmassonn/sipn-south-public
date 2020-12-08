#          Group               Forecast ID      Color for plots [string]     Comment
#          name [string]       list [list]
np.random.seed(10)
info = [ 
    [ "NASA-GSFC",        [1]               , [0.89, 0.10, 0.11], "$^{s}$"  ],  \
    [ "FIO-ESM",          [1]               , [1.00, 0.50, 0.00], "$^d$"    ],  \
#    [ "barreira",         [1]               , [0.65, 0.81, 0.89], "$^{s}$"  ],  \
    [ "NicoSun",          range(1, 3 + 1)   , [0.79, 0.70, 0.84], "$^s$"    ],  \
    [ "ucl",              range(1, 10 + 1)  , [0.42, 0.24, 0.60], "$^d$"    ],  \
#    [ "nasa-gmao",        range(1, 10 + 1)  , [0.99, 0.75, 0.44], "$^d$"    ],  \
    [ "CanCM4i", range(1, 10 + 1)  , [0.70, 0.87, 0.54], "$^d$"    ],  \
    [ "MetOffice",        range(1, 42 + 1)  , [0.98, 0.60, 0.60], "$^d$"    ],  \
    [ "CNRM",             range(1, 51 + 1)  , [0.98, 0.58, 0.88], "$^d$"    ],  \
#    [ "Lamont",           [1]               , [0.20, 0.63, 0.17], "$^{s,i}$"],  \
    [ "ecmwf",            range(1, 51 + 1)  , [0.69, 0.35, 0.16], "$^d$"    ],  \
     [ "cmcc",             range(1, 50 + 1)  , [1.00, 1.00, 0.60], "$^d$"    ],  \

       ]


