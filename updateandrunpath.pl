#!usr/bin/perl

#Written by: Wheaton Schroeder
#Written to quickly make updates to the p-ath780 model without
#running each convert file individually
system("python p-athLeaf/convertWLS.py");
system("python p-athRoot/convertWLS.py");
system("python p-athSeed/convertWLS.py");
system("python p-athStem/convertWLS.py");
system("perl makeGrowthInputs.pl");
system("gams p-ath780_revised.gms");