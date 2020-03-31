#!usr/bin/perl

#Written by: Wheaton Schroeder
#designed to format all the outputs from kathMILPgrowthCurveD.gms and from each of the 
#four Arabidopsis tissue models into the appopriate GAMS input files for KOPTIC
#additionall it creates 10 parallel code of each KOPTIC version, for a total of 90 codes to run
#each code starts at a different point in the reaction list, they are staggered by 90

use strict;

#read the pool file
open(CONCMAT, "<in_silico_conc.txt") or die "could not open the in silico concentration file, reason: $!";

#read the reaction rate file
open(RXNMAT, "<rxn_rate_out.txt") or die "could not open the reaction rate file, reason: $!";

#put reformatted file into a new file
open(CONCOUT, ">pool_matrix.txt") or die "could not write pool matrix file, reason: $!";

#put reformatted file into a new file
open(RXNOUT, ">rxn_matrix.txt") or die "could not write reaction matrix file, reason: $!";

#read both as an array
chomp(my @poolMat = <CONCMAT>);
chomp(my @rxnMat = <RXNMAT>);

#create correct formatting for pool matrix

#loop for each line in pool 2-D matrix, finish formatting
for(my $i = 0; $i <= $#poolMat; $i++) {

	#replace quote space and quote with quote period and quote
	$poolMat[$i] =~ s/\"\s\"/\"\.\"/;
	$poolMat[$i] =~ s/\"/\'/ig;
	printf CONCOUT "%s\n", $poolMat[$i];

}

#repeat for reaction matrix

#loop for each line in pool 2-D matrix, finish formatting
for(my $j = 0; $j <= $#rxnMat; $j++) {

	#replace quote space and quote with quote period and quote
	$rxnMat[$j] =~ s/\"\s\"/\"\.\"/;
	$rxnMat[$j] =~ s/\"/\'/ig;
	printf RXNOUT "%s\n", $rxnMat[$j];

}

#create a file for the total reaction type listen
#read all the files into arrays, remove all of the "/" characters
open(TYPELEAF, "<k-athLeaf/rxntype.txt") or die "could not open the leaf reaction type file, reason: $!";
chomp(my @leafType = <TYPELEAF>);
shift(@leafType);
pop(@leafType);

open(TYPEROOT, "<k-athRoot/rxntype.txt") or die "could not open the root reaction type file, reason: $!";
chomp(my @rootType = <TYPEROOT>);
shift(@rootType);
pop(@rootType);

open(TYPESEED, "<k-athSeed/rxntype.txt") or die "could not open the seed reaction type file, reason: $!";
chomp(my @seedType = <TYPESEED>);
shift(@seedType);
pop(@seedType);

open(TYPESTEM, "<k-athStem/rxntype.txt") or die "could not open the stem reaction type file, reason: $!";
chomp(my @stemType = <TYPESTEM>);
shift(@stemType);
pop(@stemType);

#write the new combined reaction file
open(TYPEALL, ">rxntype_all.txt") or die "could not write the reaction type file for all reactions, reason: $!";

printf TYPEALL "\/\n";

for(my $k = 0; $k <= $#leafType; $k++) {

	printf TYPEALL "%s\n", $leafType[$k];

}

for(my $m = 0; $m <= $#rootType; $m++) {

	printf TYPEALL "%s\n", $rootType[$m];

}

for(my $n = 0; $n <= $#seedType; $n++) {

	printf TYPEALL "%s\n", $seedType[$n];

}

for(my $p = 0; $p <= $#stemType; $p++) {

	printf TYPEALL "%s\n", $stemType[$p];

}

printf TYPEALL "\/";

printf "Reaction type files written.\n";

#create the combined Sij matrix
#start by reading all of the Sij matrix files and removing all of the "/" characters
#at the start and end of the files
open(SIJLEAF, "<k-athLeaf/S_matrix.txt") or die "could not open the leaf Sij file, reason: $!";
chomp(my @leafSij = <SIJLEAF>);
shift(@leafSij);
pop(@leafSij);

open(SIJROOT, "<k-athRoot/S_matrix.txt") or die "could not open the root Sij file, reason: $!";
chomp(my @rootSij = <SIJROOT>);
shift(@rootSij);
pop(@rootSij);

open(SIJSEED, "<k-athSeed/S_matrix.txt") or die "could not open the seed Sij file, reason: $!";
chomp(my @seedSij = <SIJSEED>);
shift(@seedSij);
pop(@seedSij);

open(SIJSTEM, "<k-athStem/S_matrix.txt") or die "could not open the stem Sij file, reason: $!";
chomp(my @stemSij = <SIJSTEM>);
shift(@stemSij);
pop(@stemSij);

my @fullSij;

#create a new Sij files for the whole system
open(SIJALL, ">S_matrix_all.txt") or die "could not write the the total Sij matrix file, reason: $!";

printf SIJALL "\/\n";

for(my $q = 0; $q <= $#leafSij; $q++) {

	printf SIJALL "%s\n", $leafSij[$q];
	push @fullSij, $leafSij[$q];

}

for(my $r = 0; $r <= $#rootSij; $r++) {

	printf SIJALL "%s\n", $rootSij[$r];
	push @fullSij, $rootSij[$r];

}

for(my $s = 0; $s <= $#seedSij; $s++) {

	printf SIJALL "%s\n", $seedSij[$s];
	push @fullSij, $seedSij[$s];

}

for(my $t = 0; $t <= $#stemSij; $t++) {

	printf SIJALL "%s\n", $stemSij[$t];
	push @fullSij, $stemSij[$t];

}

printf SIJALL "\/";

printf "S matrix files written.\n";

#create the excess metabolite files:
#excess_nr.txt -> no excess metabolites, all metabolites are allowable regulators (no restrictions)
#excess_npw.txt -> excess metabolites are protons and water (no protons or water)
#excess_npwe.txt -> excess metabolites are protons, water, and energy molecules

#start by reading the total metabolites file, remove trailing and leading "/" characters
open(ALLMETS, "<all_mets.txt") or die "could not open list of all metabolites, reason: $!";
chomp(my @allmets = <ALLMETS>);
shift(@allmets);
pop(@allmets);

#create output files
open(EXCESSNR, ">excess_nr.txt") or die "could not create excess_nr.txt, reason: $!";

open(EXCESSNPW, ">excess_npw.txt") or die "could not create excess_npw.txt, reason: $!";

open(EXCESSNPWE, ">excess_npwe.txt") or die "could not create excess_npwe.txt, reason: $!";

printf EXCESSNR "\/\n";
printf EXCESSNPW "\/\n";
printf EXCESSNPWE "\/\n";

#create hash for energy molecules
my %energy_mol = (

	"C00001" => "yes",
	"C00002" => "yes",
	"C00003" => "yes",
	"C00004" => "yes",
	"C00005" => "yes",
	"C00006" => "yes",
	"C00015" => "yes",
	"C00016" => "yes",
	"C00020" => "yes",
	"C00035" => "yes",
	"C00044" => "yes",
	"C00063" => "yes",
	"C00075" => "yes",
	"C00080" => "yes",
	"C00105" => "yes",
	"C00112" => "yes",
	"C00459" => "yes",
	"C01352" => "yes",
	"C05359" => "yes",

);

#write the files
for (my $u = 0; $u <= $#allmets; $u++) {

	#no metabolite is in excess in excess_nr.txt
	if ($allmets[$u] =~ /(begin|include)/ig) {
		
		#if a label metabolite, count as an excess so it is not used
		printf EXCESSNR "%s 1\n", $allmets[$u];
		
	} else {
	
		#if not a label metabolite, then don't count it as in excess
		printf EXCESSNR "%s 0\n", $allmets[$u];
	
	}
	
	#water and protons are in excess in excess_npw.txt
	if ($allmets[$u] =~ /(C00001|C00080)/) {
	
		printf EXCESSNPW "%s 1\n", $allmets[$u];
	
	} elsif ($allmets[$u] =~ /(begin|include)/ig) {
	
			#if a label metabolite, count as an excess so it is not used
			printf EXCESSNPW "%s 1\n", $allmets[$u];
			#for some reason this condition is ignored by the compiler in favor of its clone in the else condition
	
	} else {
	
		if($allmets[$u] =~ /(begin|include)/ig) {
	
			#if a label metabolite, count as an excess so it is not used
			printf EXCESSNPW "%s 1\n", $allmets[$u];
			#basically this is where the label metabolites are set as excess. not sure why.
	
		} else {
		
			printf EXCESSNPW "%s 0\n", $allmets[$u];
		
		}
	
	}
	
	my $key = $allmets[$u];
	
	#remove quote marks
	$key =~ s/\'//g;
	#remove compartmentalization
	$key =~ s/\[.+\]//g;
	
	if (exists $energy_mol{$key}) {
	
		#if it is an energy molecule, then count as excess for excess_npwe.txt
		printf EXCESSNPWE "%s 1\n", $allmets[$u];
	
	} elsif ($allmets[$u] =~ /(C00001|C00080)/) {
	
		#if water or proton is in excess
		printf EXCESSNPWE "%s 1\n", $allmets[$u];
	
	} elsif ($allmets[$u] =~ /(begin|include)/ig) {
	
		#if a label metabolite, count as an excess so it is not used
		printf EXCESSNPWE "%s 1\n", $allmets[$u];
	
	} else {
	
		if($allmets[$u] =~ /(begin|include)/ig) {
	
			#if a label metabolite, count as an excess so it is not used
			printf EXCESSNPWE "%s 1\n", $allmets[$u];
			#basically this is where the label metabolites are set as excess. not sure why.
	
		} else {
		
			#not proton, water, energy molecule or label, so it is not in excess, can be considered for regulatory mechanism
			printf EXCESSNPWE "%s 0\n", $allmets[$u];
		
		}
	
	}

}

printf EXCESSNR "\/";
printf EXCESSNPW "\/";
printf EXCESSNPWE "\/";

printf "Excess files written.\n";

#array for all reactions, used to write "all_rxns_X.txt" files
my @allrxns = ( );

#create metabolite files which stores which compartment each metabolite is in
open(LEAFRXNS, "<k-athLeaf/reactions.txt") or die "could not open leaf reaction list, reason: $!";
chomp(my @leafrxns = <LEAFRXNS>);
shift(@leafrxns);
pop(@leafrxns);

open(ROOTRXNS, "<k-athRoot/reactions.txt") or die "could not open root reaction list, reason: $!";
chomp(my @rootrxns = <ROOTRXNS>);
shift(@rootrxns);
pop(@rootrxns);

open(SEEDRXNS, "<k-athSeed/reactions.txt") or die "could not open seed reaction list, reason: $!";
chomp(my @seedrxns = <SEEDRXNS>);
shift(@seedrxns);
pop(@seedrxns);

open(STEMRXNS, "<k-athStem/reactions.txt") or die "could not open stem reaction list, reason: $!";
chomp(my @stemrxns = <STEMRXNS>);
shift(@stemrxns);
pop(@stemrxns);

#write to output reaction file
open(RXNSINLEAF, ">rxns_in_leaf.txt") or die "could not write file stating which reactions are in the leaf, reason: $!";

open(RXNSINROOT, ">rxns_in_root.txt") or die "could not write file stating which reactions are in the root, reason: $!";

open(RXNSINSEED, ">rxns_in_seed.txt") or die "could not write file stating which reactions are in the seed, reason: $!";

open(RXNSINSTEM, ">rxns_in_stem.txt") or die "could not write file stating which reactions are in the stem, reason: $!";

printf RXNSINLEAF "\/\n";
printf RXNSINROOT "\/\n";
printf RXNSINSEED "\/\n";
printf RXNSINSTEM "\/\n";

for (my $a = 0; $a <= $#leafrxns; $a++) {

	printf RXNSINLEAF "%s 1\n", $leafrxns[$a];
	printf RXNSINROOT "%s 0\n", $leafrxns[$a];
	printf RXNSINSEED "%s 0\n", $leafrxns[$a];
	printf RXNSINSTEM "%s 0\n", $leafrxns[$a];
	push @allrxns, $leafrxns[$a];

}

for (my $b = 0; $b <= $#rootrxns; $b++) {

	printf RXNSINLEAF "%s 0\n", $rootrxns[$b];
	printf RXNSINROOT "%s 1\n", $rootrxns[$b];
	printf RXNSINSEED "%s 0\n", $rootrxns[$b];
	printf RXNSINSTEM "%s 0\n", $rootrxns[$b];
	push @allrxns, $rootrxns[$b];

}

for (my $c = 0; $c <= $#seedrxns; $c++) {

	printf RXNSINLEAF "%s 0\n", $seedrxns[$c];
	printf RXNSINROOT "%s 0\n", $seedrxns[$c];
	printf RXNSINSEED "%s 1\n", $seedrxns[$c];
	printf RXNSINSTEM "%s 0\n", $seedrxns[$c];
	push @allrxns, $seedrxns[$c];

}

for (my $d = 0; $d <= $#stemrxns; $d++) {

	printf RXNSINLEAF "%s 0\n", $stemrxns[$d];
	printf RXNSINROOT "%s 0\n", $stemrxns[$d];
	printf RXNSINSEED "%s 0\n", $stemrxns[$d];
	printf RXNSINSTEM "%s 1\n", $stemrxns[$d];
	push @allrxns, $stemrxns[$d];

}

printf RXNSINLEAF "\/";
printf RXNSINROOT "\/";
printf RXNSINSEED "\/";
printf RXNSINSTEM "\/";

open(ALLRXNS, ">all_rxns_1.txt") or die "could not write all_rxns_1.txt file, reason: $!";
printf ALLRXNS "\/\n";

for(my $y = 0; $y <= $#allrxns; $y++) {

	printf ALLRXNS "%s\n", $allrxns[$y];

}

printf ALLRXNS "\/";

#shifts the reaction list to result in shifted starting points
for(my $ab = 0; $ab <= 90; $ab++) {
	
	#takes the first item 
	my $temp2 = shift @allrxns;
	push @allrxns, $temp2;
	
}

#write the staggered all_rxns files. Creates 90 files in total.
#also creates new KOPTIC files to parallelize the runs
#start at 2 since the first version of each file will serve as the template

#read the template files for same t files
open(SAMETNRTEMP, "<KOPTIC_nr_same_t_1.gms") or die "could not open KOPTIC_nr_same_t template, reason: $!";
chomp(my @sametnr_temp = <SAMETNRTEMP>);
my $sameTnrTemp = join "\n",@sametnr_temp;

open(SAMETNPWTEMP, "<KOPTIC_npw_same_t_1.gms") or die "could not open KOPTIC_npw_same_t template, reason: $!";
chomp(my @sametnpw_temp = <SAMETNPWTEMP>);
my $sameTnpwTemp = join "\n",@sametnpw_temp;

open(SAMETNPWETEMP, "<KOPTIC_npwe_same_t_1.gms") or die "could not open KOPTIC_npwe_same_t template, reason: $!";
chomp(my @sametnpwe_temp = <SAMETNPWETEMP>);
my $sameTnpweTemp = join "\n",@sametnpwe_temp;

#read template files for all t files
open(ALLTNRTEMP, "<KOPTIC_nr_all_t_1.gms") or die "could not open KOPTIC_nr_all_t template, reason: $!";
chomp(my @alltnr_temp = <ALLTNRTEMP>);
my $allTnrTemp = join "\n",@alltnr_temp;

open(ALLTNPWTEMP, "<KOPTIC_npw_all_t_1.gms") or die "could not open KOPTIC_npw_all_t template, reason: $!";
chomp(my @alltnpw_temp = <ALLTNPWTEMP>);
my $allTnpwTemp = join "\n",@alltnpw_temp;

open(ALLTNPWETEMP, "<KOPTIC_npwe_all_t_1.gms") or die "could not open KOPTIC_npwe_all_t template, reason: $!";
chomp(my @alltnpwe_temp = <ALLTNPWETEMP>);
my $allTnpweTemp = join "\n",@alltnpwe_temp;

#read template files for same c files
open(SAMECNRTEMP, "<KOPTIC_nr_same_c_1.gms") or die "could not open KOPTIC_nr_same_c template, reason: $!";
chomp(my @samecnr_temp = <SAMECNRTEMP>);
my $sameCnrTemp = join "\n",@samecnr_temp;

open(SAMECNPWTEMP, "<KOPTIC_npw_same_c_1.gms") or die "could not open KOPTIC_npw_same_c template, reason: $!";
chomp(my @samecnpw_temp = <SAMECNPWTEMP>);
my $sameCnpwTemp = join "\n",@samecnpw_temp;

open(SAMECNPWETEMP, "<KOPTIC_npwe_same_c_1.gms") or die "could not open KOPTIC_npwe_same_c template, reason: $!";
chomp(my @samecnpwe_temp = <SAMECNPWETEMP>);
my $sameCnpweTemp = join "\n",@samecnpwe_temp;

for(my $v = 2; $v <= 10; $v++) {

	#temporarily stores the file text
	my $temp;

	#open the files to write to and write the templates to those files with different output file names
	
	#start with same tissue KOPTIC files
	open(NRSAMET, ">KOPTIC_nr_same_t_".$v.".gms") or die "could not write new KOPTIC same tissue nr file $v, reason: $!";
	
	$temp = $sameTnrTemp;
	$temp =~ s/_1/_$v/ig;
	printf NRSAMET "%s", $temp;
	
	open(NPWSAMET, ">KOPTIC_npw_same_t_".$v.".gms") or die "could not write new KOPTIC same tissue npw file $v, reason: $!";
	
	$temp = $sameTnpwTemp;
	$temp =~ s/_1/_$v/ig;
	printf NPWSAMET "%s", $temp;
	
	open(NPWESAMET, ">KOPTIC_npwe_same_t_".$v.".gms") or die "could not write new KOPTIC same tissue npwe file $v, reason: $!";
	
	$temp = $sameTnpweTemp;
	$temp =~ s/_1/_$v/ig;
	printf NPWESAMET "%s", $temp;
	
	#all tissues KOPTIC files
	open(NRALLT, ">KOPTIC_nr_all_t_".$v.".gms") or die "could not write new KOPTIC all tissue nr file $v, reason: $!";
	
	$temp = $allTnrTemp;
	$temp =~ s/_1/_$v/ig;
	printf NRALLT "%s", $temp;
	
	open(NPWALLT, ">KOPTIC_npw_all_t_".$v.".gms") or die "could not write new KOPTIC all tissue npw file $v, reason: $!";
	
	$temp = $allTnpwTemp;
	$temp =~ s/_1/_$v/ig;
	printf NPWALLT "%s", $temp;
	
	open(NPWEALLT, ">KOPTIC_npwe_all_t_".$v.".gms") or die "could not write new KOPTIC all tissue npwe file $v, reason: $!";
	
	$temp = $allTnpweTemp;
	$temp =~ s/_1/_$v/ig;
	printf NPWEALLT "%s", $temp;
	
	#finally write the same compartment KOPTIC files
	open(NRSAMEC, ">KOPTIC_nr_same_c_".$v.".gms") or die "could not write new KOPTIC same compartment nr file $v, reason: $!";
	
	$temp = $sameCnrTemp;
	$temp =~ s/_1/_$v/ig;
	printf NRSAMEC "%s", $temp;
	
	open(NPWSAMEC, ">KOPTIC_npw_same_c_".$v.".gms") or die "could not write new KOPTIC same compartment npw file $v, reason: $!";
	
	$temp = $sameCnpwTemp;
	$temp =~ s/_1/_$v/ig;
	printf NPWSAMEC "%s", $temp;
	
	open(NPWESAMEC, ">KOPTIC_npwe_same_c_".$v.".gms") or die "could not write new KOPTIC same compartment npwe file $v, reason: $!";
	
	$temp = $sameCnpweTemp;
	$temp =~ s/_1/_$v/ig;
	printf NPWESAMEC "%s", $temp;
	
	#now write the corresponding reaction list file
	
	
	open(ALLRXNSTEMP, ">all_rxns_".$v.".txt") or die "could not write all_rxns_$v.txt file, reason: $!";
	printf ALLRXNSTEMP "\/\n";

	for(my $z = 0; $z <= $#allrxns; $z++) {

		printf ALLRXNSTEMP "%s\n", $allrxns[$z];

	}

	printf ALLRXNSTEMP "\/";
	
	#shift all_rxns by 10 to stagger starting points
	
	for(my $aa = 0; $aa <= 90; $aa++) {
	
		#takes the first item 
		my $temp3 = shift @allrxns;
		push @allrxns, $temp3;
	
	}
	
}

printf "Parallel GAMS code written.\n";
printf "Reaction list files written.\n";

#repeat the above procedure for metabolites

#read in all the metabolite files, remove starting and trailing "/" characters on the files
open(LEAFMETS, "<k-athLeaf/metabolites.txt") or die "could not open leaf metabolite list, reason: $!";
chomp(my @leafmets = <LEAFMETS>);
shift(@leafmets);
pop(@leafmets);

open(ROOTMETS, "<k-athRoot/metabolites.txt") or die "could not open root metabolite list, reason: $!";
chomp(my @rootmets = <ROOTMETS>);
shift(@rootmets);
pop(@rootmets);

open(SEEDMETS, "<k-athSeed/metabolites.txt") or die "could not open seed metabolite list, reason: $!";
chomp(my @seedmets = <SEEDMETS>);
shift(@seedmets);
pop(@seedmets);

open(STEMMETS, "<k-athStem/metabolites.txt") or die "could not open stem metabolite list, reason: $!";
chomp(my @stemmets = <STEMMETS>);
shift(@stemmets);
pop(@stemmets);

#write to output metabolite file
open(METSINLEAF, ">mets_in_leaf.txt") or die "could not write file indicating which metabolites are in the leaf, reason: $!";

open(METSINROOT, ">mets_in_root.txt") or die "could not write file indicating which metabolites are in the root, reason: $!";

open(METSINSEED, ">mets_in_seed.txt") or die "could not write file indicating which metabolites are in the seed, reason: $!";

open(METSINSTEM, ">mets_in_stem.txt") or die "could not write file indicating which metabolites are in the stem, reason: $!";

open(ALLMETS, ">all_mets.txt") or die "could not write file writing the list of all metabolites, reason: $!";

printf METSINLEAF "\/\n";
printf METSINROOT "\/\n";
printf METSINSEED "\/\n";
printf METSINSTEM "\/\n";
printf ALLMETS "\/\n";

for (my $w = 0; $w <= $#leafmets; $w++) {

	printf METSINLEAF "%s 1\n", $leafmets[$w];
	printf METSINROOT "%s 0\n", $leafmets[$w];
	printf METSINSEED "%s 0\n", $leafmets[$w];
	printf METSINSTEM "%s 0\n", $leafmets[$w];
	printf ALLMETS "%s\n", $leafmets[$w];

}

for (my $x = 0; $x <= $#rootmets; $x++) {

	printf METSINLEAF "%s 0\n", $rootmets[$x];
	printf METSINROOT "%s 1\n", $rootmets[$x];
	printf METSINSEED "%s 0\n", $rootmets[$x];
	printf METSINSTEM "%s 0\n", $rootmets[$x];
	printf ALLMETS "%s\n", $rootmets[$x];

}

for (my $y = 0; $y <= $#seedmets; $y++) {

	printf METSINLEAF "%s 0\n", $seedmets[$y];
	printf METSINROOT "%s 0\n", $seedmets[$y];
	printf METSINSEED "%s 1\n", $seedmets[$y];
	printf METSINSTEM "%s 0\n", $seedmets[$y];
	printf ALLMETS "%s\n", $seedmets[$y];

}

for (my $z = 0; $z <= $#stemmets; $z++) {

	printf METSINLEAF "%s 0\n", $stemmets[$z];
	printf METSINROOT "%s 0\n", $stemmets[$z];
	printf METSINSEED "%s 0\n", $stemmets[$z];
	printf METSINSTEM "%s 1\n", $stemmets[$z];
	printf ALLMETS "%s\n", $stemmets[$z];

}

printf METSINLEAF "\/";
printf METSINROOT "\/";
printf METSINSEED "\/";
printf METSINSTEM "\/";
printf ALLMETS "\/";


#write cell-compartment files
#read list of all reactions

printf "making reaction compartmentalization files.\n";

open(FULLRXNLIST, "<all_rxns_1.txt") or die "could not open the pool size file, reason: $!";

chomp(my @RxnList = <FULLRXNLIST>);

#remove leading and trailing "/" symbols
shift @RxnList;
pop @RxnList;

#open files to write to for specific compartments
#leaf compartments
open(RXNLT, ">rxns_in_leaf_t.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNLCL, ">rxns_in_leaf_cl.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNLE, ">rxns_in_leaf_e.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNLC, ">rxns_in_leaf_c.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNLP, ">rxns_in_leaf_p.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNLOM, ">rxns_in_leaf_om.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNLIM, ">rxns_in_leaf_im.txt") or die "could not write/create compartmentalization file, reason: $!";

#root compartments
open(RXNRE, ">rxns_in_root_e.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNRC, ">rxns_in_root_c.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNRP, ">rxns_in_root_p.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNROM, ">rxns_in_root_om.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNRIM, ">rxns_in_root_im.txt") or die "could not write/create compartmentalization file, reason: $!";

#seed compartments
open(RXNSEE, ">rxns_in_seed_e.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNSEC, ">rxns_in_seed_c.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNSEP, ">rxns_in_seed_p.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNSEOM, ">rxns_in_seed_om.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNSEIM, ">rxns_in_seed_im.txt") or die "could not write/create compartmentalization file, reason: $!";

#stem compartments
open(RXNSTE, ">rxns_in_stem_e.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNSTC, ">rxns_in_stem_c.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNSTP, ">rxns_in_stem_p.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNSTOM, ">rxns_in_stem_om.txt") or die "could not write/create compartmentalization file, reason: $!";
open(RXNSTIM, ">rxns_in_stem_im.txt") or die "could not write/create compartmentalization file, reason: $!";

#write leading "/"
printf RXNLT "\/\n";
printf RXNLCL "\/\n";
printf RXNLE "\/\n";
printf RXNLC "\/\n";
printf RXNLP "\/\n";
printf RXNLOM "\/\n";
printf RXNLIM "\/\n";
printf RXNRE "\/\n";
printf RXNRC "\/\n";
printf RXNRP "\/\n";
printf RXNROM "\/\n";
printf RXNRIM "\/\n";
printf RXNSEE "\/\n";
printf RXNSEC "\/\n";
printf RXNSEP "\/\n";
printf RXNSEOM "\/\n";
printf RXNSEIM "\/\n";
printf RXNSTE "\/\n";
printf RXNSTC "\/\n";
printf RXNSTP "\/\n";
printf RXNSTOM "\/\n";
printf RXNSTIM "\/\n";

my $matchLT = 0;
my $matchLCL = 0;
my $matchLE = 0;
my $matchLC = 0;
my $matchLP = 0;
my $matchLOM = 0;
my $matchLIM = 0;
my $matchRE = 0;
my $matchRC = 0;
my $matchRP = 0;
my $matchROM = 0;
my $matchRIM = 0;
my $matchSEE = 0;
my $matchSEC = 0;
my $matchSEP = 0;
my $matchSEOM = 0;
my $matchSEIM = 0;
my $matchSTE = 0;
my $matchSTC = 0;
my $matchSTP = 0;
my $matchSTOM = 0;
my $matchSTIM = 0;

for(my $aa = 0; $aa <= $#RxnList; $aa++) {

	#check which compartment it is in, assign a 1 to that file and a zero to all others
	#do this by searching the Sij matrix entries of the given reaction and checking what compartments are in the reaction 
	#this way we can correctly compartmentalize the transport reactions as well
	
	#value of match
	$matchLT = 0;
	$matchLCL = 0;
	$matchLE = 0;
	$matchLC = 0;
	$matchLP = 0;
	$matchLOM = 0;
	$matchLIM = 0;
	$matchRE = 0;
	$matchRC = 0;
	$matchRP = 0;
	$matchROM = 0;
	$matchRIM = 0;
	$matchSEE = 0;
	$matchSEC = 0;
	$matchSEP = 0;
	$matchSEOM = 0;
	$matchSEIM = 0;
	$matchSTE = 0;
	$matchSTC = 0;
	$matchSTP = 0;
	$matchSTOM = 0;
	$matchSTIM = 0;
	
	#format search string properly
	my $temp4 = $RxnList[$aa];
	$temp4 =~ s/\[/\\[/;
	$temp4 =~ s/\]/\\]/;
	
	#skip label reactions
	if ($RxnList[$aa] =~ /Rlabel/) {
	
		#printe the label and a zero for all compartments (don't treat labels as in a compartment)
		printf RXNLT "%s %s\n", $RxnList[$aa], $matchLT;
		printf RXNLCL "%s %s\n", $RxnList[$aa], $matchLCL;
		printf RXNLE "%s %s\n", $RxnList[$aa], $matchLE;
		printf RXNLC "%s %s\n", $RxnList[$aa], $matchLC;
		printf RXNLP "%s %s\n", $RxnList[$aa], $matchLP;
		printf RXNLOM "%s %s\n", $RxnList[$aa], $matchLOM;
		printf RXNLIM "%s %s\n", $RxnList[$aa], $matchLIM;
		printf RXNRE "%s %s\n", $RxnList[$aa], $matchRE;
		printf RXNRC "%s %s\n", $RxnList[$aa], $matchRC;
		printf RXNRP "%s %s\n", $RxnList[$aa], $matchRP;
		printf RXNROM "%s %s\n", $RxnList[$aa], $matchROM;
		printf RXNRIM "%s %s\n", $RxnList[$aa], $matchRIM;
		printf RXNSEE "%s %s\n", $RxnList[$aa], $matchSEE;
		printf RXNSEC "%s %s\n", $RxnList[$aa], $matchSEC;
		printf RXNSEP "%s %s\n", $RxnList[$aa], $matchSEP;
		printf RXNSEOM "%s %s\n", $RxnList[$aa], $matchSEOM;
		printf RXNSEIM "%s %s\n", $RxnList[$aa], $matchSEIM;
		printf RXNSTE "%s %s\n", $RxnList[$aa], $matchSTE;
		printf RXNSTC "%s %s\n", $RxnList[$aa], $matchSTC;
		printf RXNSTP "%s %s\n", $RxnList[$aa], $matchSTP;
		printf RXNSTOM "%s %s\n", $RxnList[$aa], $matchSTOM;
		printf RXNSTIM "%s %s\n", $RxnList[$aa], $matchSTIM;
	
	} else {
	
		for(my $ab = 0; $ab <= $#fullSij; $ab++) {
	
			#if the Sij line is for the current reaction
			if($fullSij[$ab] =~ /$temp4/) {
		
				#check what compartment the metabolite is include
				#if in leaf thylakoid
				if($fullSij[$ab] =~ /\[l_t\]/) {
				
					#mark that a match has been made
					$matchLT = 1;
				
				}
				
				#if in leaf chloroplast etc. 
				if($fullSij[$ab] =~ /\[l_cl\]/) {
				
					#mark that a match has been made
					$matchLCL = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[l_e\]/) {
				
					#mark that a match has been made
					$matchLE = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[l_c\]/) {
				
					#mark that a match has been made
					$matchLC = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[l_p\]/) {
				
					#mark that a match has been made
					$matchLP = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[l_om\]/) {
				
					#mark that a match has been made
					$matchLOM = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[l_im\]/) {
				
					#mark that a match has been made
					$matchLIM = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[r_e\]/) {
				
					#mark that a match has been made
					$matchRE = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[r_c\]/) {
				
					#mark that a match has been made
					$matchRC = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[r_p\]/) {
				
					#mark that a match has been made
					$matchRP = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[r_om\]/) {
				
					#mark that a match has been made
					$matchROM = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[r_im\]/) {
				
					#mark that a match has been made
					$matchRIM = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[se_e\]/) {
				
					#mark that a match has been made
					$matchSEE = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[se_c\]/) {
				
					#mark that a match has been made
					$matchSEC = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[se_p\]/) {
				
					#mark that a match has been made
					$matchSEP = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[se_om\]/) {
				
					#mark that a match has been made
					$matchSEOM = 1;
				
				}
				
				if($fullSij[$ab] =~ /\[se_im\]/) {
				
					#mark that a match has been made
					$matchSEIM = 1;
				
				}
		
			}
	
		}
		
		#prints resutls to all files, 1 if there is a match (i.e. reaction involves that compartment), 0 otherwise
		printf RXNLT "%s %s\n", $RxnList[$aa], $matchLT;
		printf RXNLCL "%s %s\n", $RxnList[$aa], $matchLCL;
		printf RXNLE "%s %s\n", $RxnList[$aa], $matchLE;
		printf RXNLC "%s %s\n", $RxnList[$aa], $matchLC;
		printf RXNLP "%s %s\n", $RxnList[$aa], $matchLP;
		printf RXNLOM "%s %s\n", $RxnList[$aa], $matchLOM;
		printf RXNLIM "%s %s\n", $RxnList[$aa], $matchLIM;
		printf RXNRE "%s %s\n", $RxnList[$aa], $matchRE;
		printf RXNRC "%s %s\n", $RxnList[$aa], $matchRC;
		printf RXNRP "%s %s\n", $RxnList[$aa], $matchRP;
		printf RXNROM "%s %s\n", $RxnList[$aa], $matchROM;
		printf RXNRIM "%s %s\n", $RxnList[$aa], $matchRIM;
		printf RXNSEE "%s %s\n", $RxnList[$aa], $matchSEE;
		printf RXNSEC "%s %s\n", $RxnList[$aa], $matchSEC;
		printf RXNSEP "%s %s\n", $RxnList[$aa], $matchSEP;
		printf RXNSEOM "%s %s\n", $RxnList[$aa], $matchSEOM;
		printf RXNSEIM "%s %s\n", $RxnList[$aa], $matchSEIM;
		printf RXNSTE "%s %s\n", $RxnList[$aa], $matchSTE;
		printf RXNSTC "%s %s\n", $RxnList[$aa], $matchSTC;
		printf RXNSTP "%s %s\n", $RxnList[$aa], $matchSTP;
		printf RXNSTOM "%s %s\n", $RxnList[$aa], $matchSTOM;
		printf RXNSTIM "%s %s\n", $RxnList[$aa], $matchSTIM;
		
	}

}

#write trailing "/"
printf RXNLT "\/";
printf RXNLCL "\/";
printf RXNLE "\/";
printf RXNLC "\/";
printf RXNLP "\/";
printf RXNLOM "\/";
printf RXNLIM "\/";
printf RXNRE "\/";
printf RXNRC "\/";
printf RXNRP "\/";
printf RXNROM "\/";
printf RXNRIM "\/";
printf RXNSEE "\/";
printf RXNSEC "\/";
printf RXNSEP "\/";
printf RXNSEOM "\/";
printf RXNSEIM "\/";
printf RXNSTE "\/";
printf RXNSTC "\/";
printf RXNSTP "\/";
printf RXNSTOM "\/";
printf RXNSTIM "\/";

printf "making metabolite compartmentalization files.\n";

open(FULLMETLIST, "<all_mets.txt") or die "could not open the list of all metabolites, reason: $!";

chomp(my @MetList = <FULLMETLIST>);

#remove leading and trailing "/" symbols
shift @MetList;
pop @MetList;

#open files to write to for specific compartments
#leaf compartments
open(METLT, ">mets_in_leaf_t.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METLCL, ">mets_in_leaf_cl.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METLE, ">mets_in_leaf_e.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METLC, ">mets_in_leaf_c.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METLP, ">mets_in_leaf_p.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METLOM, ">mets_in_leaf_om.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METLIM, ">mets_in_leaf_im.txt") or die "could not write/create compartmentalization file, reason: $!";

#root compartments
open(METRE, ">mets_in_root_e.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METRC, ">mets_in_root_c.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METRP, ">mets_in_root_p.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METROM, ">mets_in_root_om.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METRIM, ">mets_in_root_im.txt") or die "could not write/create compartmentalization file, reason: $!";

#seed compartments
open(METSEE, ">mets_in_seed_e.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METSEC, ">mets_in_seed_c.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METSEP, ">mets_in_seed_p.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METSEOM, ">mets_in_seed_om.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METSEIM, ">mets_in_seed_im.txt") or die "could not write/create compartmentalization file, reason: $!";

#stem compartments
open(METSTE, ">mets_in_stem_e.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METSTC, ">mets_in_stem_c.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METSTP, ">mets_in_stem_p.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METSTOM, ">mets_in_stem_om.txt") or die "could not write/create compartmentalization file, reason: $!";
open(METSTIM, ">mets_in_stem_im.txt") or die "could not write/create compartmentalization file, reason: $!";

#write leading "/"
printf METLT "\/\n";
printf METLCL "\/\n";
printf METLE "\/\n";
printf METLC "\/\n";
printf METLP "\/\n";
printf METLOM "\/\n";
printf METLIM "\/\n";
printf METRE "\/\n";
printf METRC "\/\n";
printf METRP "\/\n";
printf METROM "\/\n";
printf METRIM "\/\n";
printf METSEE "\/\n";
printf METSEC "\/\n";
printf METSEP "\/\n";
printf METSEOM "\/\n";
printf METSEIM "\/\n";
printf METSTE "\/\n";
printf METSTC "\/\n";
printf METSTP "\/\n";
printf METSTOM "\/\n";
printf METSTIM "\/\n";

for(my $ac = 0; $ac <= $#MetList; $ac++) {

	#check which compartment it is in, assign a 1 to that file and a zero to all others
	#do this by searching the Sij matrix entries of the given reaction and checking what compartments are in the reaction 
	#this way we can correctly compartmentalize the transport reactions as well
	
	#value of match
	$matchLT = 0;
	$matchLCL = 0;
	$matchLE = 0;
	$matchLC = 0;
	$matchLP = 0;
	$matchLOM = 0;
	$matchLIM = 0;
	$matchRE = 0;
	$matchRC = 0;
	$matchRP = 0;
	$matchROM = 0;
	$matchRIM = 0;
	$matchSEE = 0;
	$matchSEC = 0;
	$matchSEP = 0;
	$matchSEOM = 0;
	$matchSEIM = 0;
	$matchSTE = 0;
	$matchSTC = 0;
	$matchSTP = 0;
	$matchSTOM = 0;
	$matchSTIM = 0;
	
	#format search string properly
	my $temp4 = $RxnList[$ac];
	$temp4 =~ s/\[/\\[/;
	$temp4 =~ s/\]/\\]/;
	
	#skip label reactions
	if ($MetList[$ac] =~ /include/) {
	
		#printe the label and a zero for all compartments (don't treat labels as in a compartment)
		printf METLT "%s %s\n", $MetList[$ac], $matchLT;
		printf METLCL "%s %s\n", $MetList[$ac], $matchLCL;
		printf METLE "%s %s\n", $MetList[$ac], $matchLE;
		printf METLC "%s %s\n", $MetList[$ac], $matchLC;
		printf METLP "%s %s\n", $MetList[$ac], $matchLP;
		printf METLOM "%s %s\n", $MetList[$ac], $matchLOM;
		printf METLIM "%s %s\n", $MetList[$ac], $matchLIM;
		printf METRE "%s %s\n", $MetList[$ac], $matchRE;
		printf METRC "%s %s\n", $MetList[$ac], $matchRC;
		printf METRP "%s %s\n", $MetList[$ac], $matchRP;
		printf METROM "%s %s\n", $MetList[$ac], $matchROM;
		printf METRIM "%s %s\n", $MetList[$ac], $matchRIM;
		printf METSEE "%s %s\n", $MetList[$ac], $matchSEE;
		printf METSEC "%s %s\n", $MetList[$ac], $matchSEC;
		printf METSEP "%s %s\n", $MetList[$ac], $matchSEP;
		printf METSEOM "%s %s\n", $MetList[$ac], $matchSEOM;
		printf METSEIM "%s %s\n", $MetList[$ac], $matchSEIM;
		printf METSTE "%s %s\n", $MetList[$ac], $matchSTE;
		printf METSTC "%s %s\n", $MetList[$ac], $matchSTC;
		printf METSTP "%s %s\n", $MetList[$ac], $matchSTP;
		printf METSTOM "%s %s\n", $MetList[$ac], $matchSTOM;
		printf METSTIM "%s %s\n", $MetList[$ac], $matchSTIM;
	
	} else {
		
		#check what compartment the metabolite is include
		#if in leaf thylakoid
		if($MetList[$ac] =~ /\[l_t\]/) {
				
			#mark that a match has been made
			$matchLT = 1;
				
		}
				
		#if in leaf chloroplast etc. 
		if($MetList[$ac] =~ /\[l_cl\]/) {
				
			#mark that a match has been made
			$matchLCL = 1;
				
		}
				
		if($MetList[$ac] =~ /\[l_e\]/) {
				
			#mark that a match has been made
			$matchLE = 1;
				
		}
				
		if($MetList[$ac] =~ /\[l_c\]/) {
				
			#mark that a match has been made
			$matchLC = 1;
				
		}
				
		if($MetList[$ac] =~ /\[l_p\]/) {
				
			#mark that a match has been made
			$matchLP = 1;
				
		}
				
		if($MetList[$ac] =~ /\[l_om\]/) {
				
			#mark that a match has been made
			$matchLOM = 1;
				
		}
				
		if($MetList[$ac] =~ /\[l_im\]/) {
				
			#mark that a match has been made
			$matchLIM = 1;
				
		}
				
		if($MetList[$ac] =~ /\[r_e\]/) {
				
			#mark that a match has been made
			$matchRE = 1;
				
		}
				
		if($MetList[$ac] =~ /\[r_c\]/) {
				
			#mark that a match has been made
			$matchRC = 1;
				
		}
			
		if($MetList[$ac] =~ /\[r_p\]/) {
				
			#mark that a match has been made
			$matchRP = 1;
			
		}
				
		if($MetList[$ac] =~ /\[r_om\]/) {
				
			#mark that a match has been made
			$matchROM = 1;
				
		}
				
		if($MetList[$ac] =~ /\[r_im\]/) {
				
			#mark that a match has been made
			$matchRIM = 1;
				
		}
				
		if($MetList[$ac] =~ /\[se_e\]/) {
			
			#mark that a match has been made
			$matchSEE = 1;
				
		}
				
		if($MetList[$ac] =~ /\[se_c\]/) {
			
			#mark that a match has been made
			$matchSEC = 1;
				
		}
				
		if($MetList[$ac] =~ /\[se_p\]/) {
				
			#mark that a match has been made
			$matchSEP = 1;
		
		}
				
		if($MetList[$ac] =~ /\[se_om\]/) {
				
			#mark that a match has been made
			$matchSEOM = 1;
				
		}
				
		if($MetList[$ac] =~ /\[se_im\]/) {
			
			#mark that a match has been made
			$matchSEIM = 1;
				
		}
		
		#prints resutls to all files, 1 if there is a match (i.e. reaction involves that compartment), 0 otherwise
		printf METLT "%s %s\n", $MetList[$ac], $matchLT;
		printf METLCL "%s %s\n", $MetList[$ac], $matchLCL;
		printf METLE "%s %s\n", $MetList[$ac], $matchLE;
		printf METLC "%s %s\n", $MetList[$ac], $matchLC;
		printf METLP "%s %s\n", $MetList[$ac], $matchLP;
		printf METLOM "%s %s\n", $MetList[$ac], $matchLOM;
		printf METLIM "%s %s\n", $MetList[$ac], $matchLIM;
		printf METRE "%s %s\n", $MetList[$ac], $matchRE;
		printf METRC "%s %s\n", $MetList[$ac], $matchRC;
		printf METRP "%s %s\n", $MetList[$ac], $matchRP;
		printf METROM "%s %s\n", $MetList[$ac], $matchROM;
		printf METRIM "%s %s\n", $MetList[$ac], $matchRIM;
		printf METSEE "%s %s\n", $MetList[$ac], $matchSEE;
		printf METSEC "%s %s\n", $MetList[$ac], $matchSEC;
		printf METSEP "%s %s\n", $MetList[$ac], $matchSEP;
		printf METSEOM "%s %s\n", $MetList[$ac], $matchSEOM;
		printf METSEIM "%s %s\n", $MetList[$ac], $matchSEIM;
		printf METSTE "%s %s\n", $MetList[$ac], $matchSTE;
		printf METSTC "%s %s\n", $MetList[$ac], $matchSTC;
		printf METSTP "%s %s\n", $MetList[$ac], $matchSTP;
		printf METSTOM "%s %s\n", $MetList[$ac], $matchSTOM;
		printf METSTIM "%s %s\n", $MetList[$ac], $matchSTIM;
		
	}

}

#write trailing "/"
printf RXNLT "\/";
printf RXNLCL "\/";
printf RXNLE "\/";
printf RXNLC "\/";
printf RXNLP "\/";
printf RXNLOM "\/";
printf RXNLIM "\/";
printf RXNRE "\/";
printf RXNRC "\/";
printf RXNRP "\/";
printf RXNROM "\/";
printf RXNRIM "\/";
printf RXNSEE "\/";
printf RXNSEC "\/";
printf RXNSEP "\/";
printf RXNSEOM "\/";
printf RXNSEIM "\/";
printf RXNSTE "\/";
printf RXNSTC "\/";
printf RXNSTP "\/";
printf RXNSTOM "\/";
printf RXNSTIM "\/";

printf "Files written.\nSUCCESS\n";