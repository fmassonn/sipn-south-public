#          Group     Nb forecasts  Color for plots     Comment
#          name
np.random.seed(10)
info = [ ["nrl",        [t for t in range(1, 10 + 1) if t!= 5],  "#92B558",   ""              ],\
         ["Lamont",     [1],                                     "#9C9A40",   "(interpolated)"],\
         ["nasa-gmao",  range(1, 10 + 1),                        "#4F84C4",   ""              ],\
         ["MetOffice",  range(1, 42 + 1),                        "#D2691E",   ""              ],\
         ["Nico-Sun",   range(1, 3 + 1),                         "#DC4C46",   ""              ],\
       ]



