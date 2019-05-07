#          Group               Forecast ID      Color for plots [string]     Comment
#          name [string]       list [list]
np.random.seed(10)
info = [ 
         [ "Petty-NASA",       [1]               , [0.89, 0.42, 0.04],  "$^{s,i}$"  ],  \
         [ "FIO-ESM",          [1]               , [0.89, 0.42, 0.04],  "$^d$"                ],  \
         [ "Lamont",           [1]               , [0.89, 0.42, 0.04],  "$^{s,i}$"  ],  \
         [ "barreira",         [1]               , [0.89, 0.42, 0.04],  "$^{s,i}$"  ],  \
         [ "Nico-Sun",         range(1, 3 + 1)   , [0.89, 0.42, 0.04],  "$^s$"                ],  \
         [ "ucl",              range(1, 10 + 1)  , [0.00, 0.44, 0.75],  "$^d$"                ],  \
         [ "nasa-gmao",        range(1, 10 + 1)  , [0.00, 0.44, 0.75],  "$^d$"                ],  \
         [ "nrl",              range(1, 10 + 1)  , [0.00, 0.44, 0.75],  "$^d$"                ],  \
         [ "Modified-CanSIPS", range(1, 20 + 1)  , [0.00, 0.44, 0.75],    "$^d$"                ],  \
         [ "MetOffice",        range(1, 42 + 1)  , [0.00, 0.44, 0.75],  "$^d$"                ],  \
         [ "ecmwf",            range(1, 50 + 1)  , [0.00, 0.44, 0.75],   "$^d$"                ],  \
         [ "CMCC",             range(1, 50 + 1)  , [0.00, 0.44, 0.75],  "$^d$"                ],  \
       ]



