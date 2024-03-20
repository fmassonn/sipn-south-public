#          Group               Forecast ID      Color for plots [string]     Comment
#          name [string]       list [list]
np.random.seed(10)
info = [ 
    [ "climatology",      range(1, 30  + 1),  ""  ],
    [ "AWI-CPS",          range(1, 30  + 1), "$^d",    ],
    [ "BSC",              range(1, 30 + 1),  "$^d$"    ], \
    [ "IOCAS-SIPNet",     [1],               "$^s"     ], \
    [ "Lamont",           [1],               "$^{s,i}$"], \
    [ "MetOffice",        range(1, 42 + 1),  "$^d$"    ], \
    [ "NicoSun",          range(1, 3 + 1),   "$^s$"    ], \
    [ "SYSU-SML-ConvLSTM",[1],               "$^s$"    ], \
    [ "SYSU-SML-MLM",     [1],               "$^s$"    ], \
    [ "cmcc",             range(1, 50 + 1),  "$^d$"    ], \
    [ "emc",              range(1, 10 + 1),  "$^d$"    ], \
    [ "gfdl",             range(1, 30 + 1),  "$^d$"    ], \
    [ "ucl",              range(1, 10 + 1),  "$^d$"    ], \
       ]


