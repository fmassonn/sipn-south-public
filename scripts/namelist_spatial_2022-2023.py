#          Group     Nb forecasts  Color for plots     Comment
#          name
np.random.seed(10)
info = [ 
          ["MetOffice",  range(1, 42 + 1), "#0077D4",   "$^d$"    ],  \
#          ["ucl",        range(1, 10 + 1), "#005BC3",   "$^d$"    ],  \
#          ["CNRM",       range(1, 51 + 1), "#00A2E1",   "$^d$"    ],  \
          ["NicoSun",    range(1, 3 + 1),  "#CA4E00",   "$^s$"    ],  \
#          ["barreira",   [1]             , "#E76F00",   "$^s$"    ],  \
#          ["Lamont",     [1],              "#FF8642",   "$^{s,i}$"],  \
          ["cmcc"  ,     range(1, 50 + 1), "#0070B2",   "$^d$"    ],  \
          [ "gfdl"   ,   range(1, 30 + 1), "#0047B6",   "$^d$"    ],  \
#          ["SYSU"      , [1]             , "#FF963B",   "$^s$"    ],  \
           ["SYSU_SML_ConvLSTM", [1],              "#12AD2B", "$^s$"    ],
       ]
