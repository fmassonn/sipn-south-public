#          Group     Nb forecasts  Color for plots     Comment
#          name
np.random.seed(10)
info = [ 
         ["climatology", range(1, 30 + 1), "grey", ""], \
         ["nasa-gmao",  range(1, 10 + 1), [0.00, 0.44, 0.75],   "$^d$"   ],  \
         ["MetOffice",  range(1, 42 + 1), [0.00, 0.64, 0.85],   "$^d$"   ],  \
         ["ucl",        range(1, 10 + 1), [0.00, 0.84, 0.95],   "$^d$"   ],  \
         ["CNRM",       range(1, 51 + 1), [0.00, 0.94, 0.95],   "$^d$"   ],  \
         ["NicoSun",    range(1, 3 + 1),  [0.89, 0.42, 0.04],   "$^s$"   ],  \
         ["barreira",   [1]             , [0.89, 0.62, 0.24],   "$^s$"   ],  \
         ["Lamont",     [1],              [0.89, 0.82, 0.44],   "$^{s,i}$"], \

       ]

