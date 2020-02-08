#          Group               Forecast ID      Color for plots [string]     Comment
#          name [string]       list [list]
np.random.seed(10)
info = [ 
         [ "NASA-GSFC",       [1]               , [0.8901960784313725, 0.10196078431372549, 0.10980392156862745],  "$^{s}$"  ],  \
         [ "FIO-ESM",          [1]               , [1.0, 0.49803921568627452, 0.0],  "$^d$"                ],  \
#         [ "Lamont",           [1]               , [0.20000000000000001, 0.62745098039215685, 0.17254901960784313],  "$^{s,i}$"  ],  \
         [ "barreira",         [1]               , [0.65098039215686276, 0.80784313725490198, 0.8901960784313725],  "$^{s}$"  ],  \
         [ "NicoSun",         range(1, 3 + 1)   , [0.792156862745098, 0.69803921568627447, 0.83921568627450982],   "$^s$"                ],  \
         [ "ucl",              range(1, 10 + 1)  , [0.41568627450980394, 0.23921568627450981, 0.60392156862745094],  "$^d$"                ],  \
         [ "nasa-gmao",        range(1, 10 + 1)  , [0.99215686274509807, 0.74901960784313726, 0.43529411764705883],  "$^d$"                ],  \
#         [ "nrl",              range(1, 10 + 1)  , [0.12156862745098039, 0.47058823529411764, 0.70588235294117652],  "$^d$"                ],  \
#         [ "Modified-CanSIPS", range(1, 20 + 1)  , [0.69803921568627447, 0.87450980392156863, 0.54117647058823526],    "$^d$"                ],  \
         [ "MetOffice",        range(1, 42 + 1)  , [0.98431372549019602, 0.60392156862745094, 0.59999999999999998],  "$^d$"                ],  \
#         [ "ecmwf",            range(1, 50 + 1)  , [0.69411764705882351, 0.34901960784313724, 0.15686274509803921],   "$^d$"                ],  \
#         [ "CMCC",             range(1, 50 + 1)  , [1.0, 1.0, 0.59999999999999998],  "$^d$"                ],  \
       ]


