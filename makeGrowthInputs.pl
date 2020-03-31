#!usr/bin/perl

#Written by: Wheaton Schroeder
#Written to create input files for p-ath780 model so these don't have to be 
#made by hand. Just making input files that change with the models, where manual
#input is not required. Saves a lot of manual labor in the long run

#NOTES:
#this file is usually named "makeGrowthInputs.pl" in the normal workflow
#This file does first rely on the "convertWLS.py" script being run on each tissue model
#This file assumes that each set of outputs from "convertWLS.py" is contained in its own
#	tissue-specific folder named "p-athTissue"
#Running this script is a prerequisite for running "p-ath780.gms"

use strict;

#read in all the reaction files, remove starting and trailing "/" characters on the files
open(LEAFRXNS, "<p-athLeaf/reactions.txt") or die "could not open leaf reaction list, reason: $!";
chomp(my @leafrxns = <LEAFRXNS>);
shift(@leafrxns);
pop(@leafrxns);

open(ROOTRXNS, "<p-athRoot/reactions.txt") or die "could not open root reaction list, reason: $!";
chomp(my @rootrxns = <ROOTRXNS>);
shift(@rootrxns);
pop(@rootrxns);

open(SEEDRXNS, "<p-athSeed/reactions.txt") or die "could not open seed reaction list, reason: $!";
chomp(my @seedrxns = <SEEDRXNS>);
shift(@seedrxns);
pop(@seedrxns);

open(STEMRXNS, "<p-athStem/reactions.txt") or die "could not open stem reaction list, reason: $!";
chomp(my @stemrxns = <STEMRXNS>);
shift(@stemrxns);
pop(@stemrxns);

#write to output reaction file
open(ALLRXNS, ">all_rxns.txt") or die "could not write list of all reactions, reason: $!";

printf ALLRXNS "\/\n";

for (my $a = 0; $a <= $#leafrxns; $a++) {

	printf ALLRXNS "%s\n", $leafrxns[$a];

}

for (my $b = 0; $b <= $#rootrxns; $b++) {

	printf ALLRXNS "%s\n", $rootrxns[$b];

}

for (my $d = 0; $d <= $#stemrxns; $d++) {

	printf ALLRXNS "%s\n", $stemrxns[$d];

}

for (my $c = 0; $c <= $#seedrxns; $c++) {

	printf ALLRXNS "%s\n", $seedrxns[$c];

}

printf ALLRXNS "\/";

#repeat the above procedure for metabolites

#read in all the metabolite files, remove starting and trailing "/" characters on the files
open(LEAFMETS, "<p-athLeaf/metabolites.txt") or die "could not open leaf metabolite list, reason: $!";
chomp(my @leafmets = <LEAFMETS>);
shift(@leafmets);
pop(@leafmets);

open(ROOTMETS, "<p-athRoot/metabolites.txt") or die "could not open root metabolite list, reason: $!";
chomp(my @rootmets = <ROOTMETS>);
shift(@rootmets);
pop(@rootmets);

open(SEEDMETS, "<p-athSeed/metabolites.txt") or die "could not open seed metabolite list, reason: $!";
chomp(my @seedmets = <SEEDMETS>);
shift(@seedmets);
pop(@seedmets);

open(STEMMETS, "<p-athStem/metabolites.txt") or die "could not open stem metabolite list, reason: $!";
chomp(my @stemmets = <STEMMETS>);
shift(@stemmets);
pop(@stemmets);

#write to output metabolite file
open(ALLMETS, ">all_mets.txt") or die "could not write list of all metabolites, reason: $!";

printf ALLMETS "\/\n";

for (my $e = 0; $e <= $#leafmets; $e++) {

	printf ALLMETS "%s\n", $leafmets[$e];

}

for (my $f = 0; $f <= $#rootmets; $f++) {

	printf ALLMETS "%s\n", $rootmets[$f];

}

for (my $g = 0; $g <= $#seedmets; $g++) {

	printf ALLMETS "%s\n", $seedmets[$g];

}

for (my $h = 0; $h <= $#stemmets; $h++) {

	printf ALLMETS "%s\n", $stemmets[$h];

}

printf ALLMETS "\/";

#repeat the above procedure for the Sij file

#read in all the metabolite files, remove starting and trailing "/" characters on the files
open(LEAFSIJ, "<p-athLeaf/S_matrix.txt") or die "could not open leaf stoichiometric list, reason: $!";
chomp(my @leafsij = <LEAFSIJ>);
shift(@leafsij);
pop(@leafsij);

open(ROOTSIJ, "<p-athRoot/S_matrix.txt") or die "could not open root stoichiometric list, reason: $!";
chomp(my @rootsij = <ROOTSIJ>);
shift(@rootsij);
pop(@rootsij);

open(SEEDSIJ, "<p-athSeed/S_matrix.txt") or die "could not open seed stoichiometric list, reason: $!";
chomp(my @seedsij = <SEEDSIJ>);
shift(@seedsij);
pop(@seedsij);

open(STEMSIJ, "<p-athStem/S_matrix.txt") or die "could not open stem stoichiometric list, reason: $!";
chomp(my @stemsij = <STEMSIJ>);
shift(@stemsij);
pop(@stemsij);

#write to output metabolite file
open(ALLSIJ, ">all_S_matrix.txt") or die "could not write list of all stoichiometric, reason: $!";

printf ALLSIJ "\/\n";

for (my $i = 0; $i <= $#leafsij; $i++) {

	printf ALLSIJ "%s\n", $leafsij[$i];

}

for (my $j = 0; $j <= $#rootsij; $j++) {

	printf ALLSIJ "%s\n", $rootsij[$j];

}

for (my $k = 0; $k <= $#seedsij; $k++) {

	printf ALLSIJ "%s\n", $seedsij[$k];

}

for (my $m = 0; $m <= $#stemsij; $m++) {

	printf ALLSIJ "%s\n", $stemsij[$m];

}

printf ALLSIJ "\/";

#repeat the above procedure for the reaction type file

#read in all the metabolite files, remove starting and trailing "/" characters on the files
open(LEAFTYPE, "<p-athLeaf/rxntype.txt") or die "could not open leaf reaction type list, reason: $!";
chomp(my @leaftype = <LEAFTYPE>);
shift(@leaftype);
pop(@leaftype);

open(ROOTTYPE, "<p-athRoot/rxntype.txt") or die "could not open root reaction type list, reason: $!";
chomp(my @roottype = <ROOTTYPE>);
shift(@roottype);
pop(@roottype);

open(SEEDTYPE, "<p-athSeed/rxntype.txt") or die "could not open seed reaction type list, reason: $!";
chomp(my @seedtype = <SEEDTYPE>);
shift(@seedtype);
pop(@seedtype);

open(STEMTYPE, "<p-athStem/rxntype.txt") or die "could not open stem reaction type list, reason: $!";
chomp(my @stemtype = <STEMTYPE>);
shift(@stemtype);
pop(@stemtype);

#write to output metabolite file
open(ALLTYPE, ">all_rxn_types.txt") or die "could not write list of all reaction type, reason: $!";

printf ALLTYPE "\/\n";

for (my $n = 0; $n <= $#leaftype; $n++) {

	printf ALLTYPE "%s\n", $leaftype[$n];

}

for (my $p = 0; $p <= $#roottype; $p++) {

	printf ALLTYPE "%s\n", $roottype[$p];

}

for (my $q = 0; $q <= $#seedtype; $q++) {

	printf ALLTYPE "%s\n", $seedtype[$q];

}

for (my $r = 0; $r <= $#stemtype; $r++) {

	printf ALLTYPE "%s\n", $stemtype[$r];

}

printf ALLTYPE "\/";

#repeat the above procedure for the reaction type file

#read in all the metabolite files, remove starting and trailing "/" characters on the files
open(LEAFREG, "<p-athLeaf/regulated_reactions.txt") or die "could not open regulated reactions list, reason: $!";
chomp(my @leafreg = <LEAFREG>);
shift(@leafreg);
pop(@leafreg);

open(ROOTREG, "<p-athRoot/regulated_reactions.txt") or die "could not open regulated reactions list, reason: $!";
chomp(my @rootreg = <ROOTREG>);
shift(@rootreg);
pop(@rootreg);

open(SEEDREG, "<p-athSeed/regulated_reactions.txt") or die "could not open regulated reactions list, reason: $!";
chomp(my @seedreg = <SEEDREG>);
shift(@seedreg);
pop(@seedreg);

open(STEMREG, "<p-athStem/regulated_reactions.txt") or die "could not open regulated reactions list, reason: $!";
chomp(my @stemreg = <STEMREG>);
shift(@stemreg);
pop(@stemreg);

#write to output metabolite file
open(ALLREG, ">regulated_reactions.txt") or die "could not write list of all reaction type, reason: $!";

printf ALLREG "\/\n";

for (my $s = 0; $s <= $#leafreg; $s++) {

	printf ALLREG "%s\n", $leafreg[$s];

}

for (my $t = 0; $t <= $#rootreg; $t++) {

	printf ALLREG "%s\n", $rootreg[$t];

}

for (my $u = 0; $u <= $#seedreg; $u++) {

	printf ALLREG "%s\n", $seedreg[$u];

}

for (my $v = 0; $v <= $#stemreg; $v++) {

	printf ALLREG "%s\n", $stemreg[$v];

}

printf ALLREG "\/";

printf "success!\n";
