*****************************************************************************
*****************************************************************************
*****************************************************************************
*Arabidopsis multi-tissue ORKA DFBA model									*
*By: Wheaton Schroeder														*
*Research Group: PI: Rajib Saha SSBIO University of Nebraska-Lincoln		*
*Latest Version: 02/18/2020													*
*																			*
*This version of the code has been written first on 02/04/2020 when it was 	*
*rightly pointed out by paper reviewers that we were indeed using something	*
*very much like a DFBA. In response to their comments, we formalized and 	*
*named this new approach ORKA. Here, we have cleaned up the code and made it*
*easier to follow so that hopefully others might use the ORKA approach		*
*themselves, or even use this model as a basis for other plant models.		*
*****************************************************************************
*****************************************************************************
*****************************************************************************

$INLINECOM /*  */
$onlisting
$ondigit
$oneolcom
$onempty

*options for this solution file
*Options function as follows:
*ensures that eight decimal places are saved for each solution
*ensures that linear problems are solved using the CPLEX solver

OPTIONS

	decimals = 8
	lp = cplex

;

*****************************************************************************
**************************** Symbols Initialized ****************************
*****************************************************************************

SETS

	I										set of metabolites in the model
$include "all_mets.txt"

	I_prime(I)								set of metabolites which must balance
$include "must_bal.txt"

	J										set of reactions in the model
$include "all_rxns.txt"

	jseed(J)								set of reactions in the seed
$include "p-athSeed/reactions.txt"

	reg(j)									regulated reactions in the leaf
$include "regulated_reactions.txt"

	specNames								names of input specification in order
$include "growthSpecNames.txt"

	hour									time points (hours) at which to test plant model and get relevant data within the day
$include "hours.txt"

	Day										set of days in which the plant grows
$include "days.txt"

	timeData								data to save for each hour
$include "timeData.txt"

	N										set of steps in the chosen Runge-Kutta method
$include "N.txt"

	theta									set of tissue in the model
$include "theta.txt"

	Nm1(n)									set n minus the last element necessary to do the summation for Runge-Kuta method

	seed_stores(j)							set of reactions representing the seed storage
$include "seed_stores.txt"

;

*define some duplicate sets used in loops
alias(J1,J);
alias(I1,I);
alias(N1,N);

*define variable sets as "NO" for now
Nm1(N) = NO;

PARAMETERS

	/*INPUT PARAMETERS*/
	
	rxntype(j)						Reaction type for rxns of leaf tissue
$include "all_rxn_types.txt"

	S(I,J)							stoichiometric matrix for all models
$include "all_S_matrix.txt"

	specs(specNames)				specifications for the optimization problem. First item is level of seeding compared to max number of seeds (scale 0 to 1) second item is current mass of plant 
$include "growthSpecs.txt"

	sunrise(day)					sunrise time table
$include "sunrise.txt"

	sunset(day)						sunset time table
$include "sunset.txt"

	timeOfDay(hour)					time of the day using military (24 hour) time
$include "timeofday.txt"

	day_num(day)					turns day set element to a number which can be used to perform operations
$include "day_number.txt"

	a(N,N)							parameter representing the interior weight values of a Runge-Kutta method
$include "a.txt"

	b(N)							parameter representing the steps in the independent variable for a Runge-Kutta method
$include "b.txt"

	c(N)							parameter representing the weights of derivative estimates for the final estimate in a Runge-Kutta method
$include "c.txt"

	up_rate(seed_stores)			uptake limits for metabolites in the seed stores for the seedlings to use up
$include "seed_up_rate.txt"

	/*the parameters below are specified in the specs file rather than having their own file*/
	c_theta(theta)		coefficient representing the change in tissue biomass as seed mass is added

	x_leafinit			mass fraction of leaf tissues in the plant at hour 0
	
	x_rootinit			mass fraction of root tissues in the plant at hour 0
	
	x_seedinit			mass fraction of seed tissues in the plant at hour 0
	
	x_steminit			mass fraction of stem tissues in the plant at hour 0
	
	stage_0_7			time (in hours) at which stage 0.7 starts important because leaves not exposed until then
	
	stage_1_2			time (in hours) at which stage 0.7 starts important because still have seed stores until then
	
	Flowering_Start		time to first flower (seed) in hours
	
	Flowering_End		time to first flower (seed) in hours
	
	Seed_Loss_Start		time to starting to loose seeds
	
	Seed_Loss_End		time at which all seeds have been lost
	
	Max_Time			maximum time allowed for the simulation
	
	inf_const			inflation constant allows us not to worry much about 
	
	light_up			light uptake limit
	
	A_l					Amplitude of sine wave fit of diurnal starch storage pattern in leaf
	
	f_l					frequency of sine wave fit of diurnal starch storage pattern in leaf
	
	b_l					x-intercept of sine wave fit of diurnal starch storage pattern in leaf
	
	A_st1				Amplitude of sine wave fit of diurnal starch storage pattern in stem
	
	f_st1				frequency of sine wave fit of diurnal starch storage pattern in stem
	
	b_st1				x-intercept of sine wave fit of diurnal starch storage pattern in stem
	
	A_st2				Amplitude of sine wave fit of diurnal sucrose storage pattern in stem
	
	f_st2				frequency of sine wave fit of diurnal sucrose storage pattern in stem
	
	b_st2				x-intercept of sine wave fit of diurnal sucrose storage pattern in stem

	k_m					rate of biomass sink into maintenance respiration
	
	k_m_dying			rate of biomass sink into maintenance respiration when the plant is dying
	
	k_s					rate of biomass sink into senescence
	
	k_s_dying			rate of biomass sink into senescence when the plant is dying
	
	Y_G					rate of biomass produced that goes to new plant growth
	
	Seeding_Gain		rate at which seeding level grows
	
	Seeding_Loss		rate at which seeding levels decrease
	
	tau_leaf			rate of transpiration from the leafd
	
	CO2_upmax			maximum rate of carbon dioxide uptake by leaf
	
	O2_upmax			maximum rate of oxygen uptake by leaf
	
	/*FIXED PARAMETERS*/
	
	M				a very large number

	epsilon			a very small number
	
	/*PARAMETERS UPDATED EACH TIMESTEP*/
	
	/*parameters constraining reaction rates*/
	UB(J)			upper bounds on reaction rate for model reactions
	
	LB(J)			lower bounds on reaction rate for model reactions

	/*parameters which store the initial conditions*/
	x_t(theta)			mass fraction of tissue theta in the plant at current hour
	
	x_n(n,theta)		mass fraction of tissue theta at runge kutta step n
	
	x_deltat(theta)		mass fraction of tissue theta in the plant at the next hour
	
	t_0					initial time
	
	t_Deltat			time at the next time point
	
	Z_0(I)				initial metabolite concentration for the current step
	
	Z_Deltat(I)			metabolite concentration for the next step
	
	Y_0					initial whole-plant mass
	
	dY_dt				time step estimate from the Runge-Kutta method
	
	Y_Deltat			whole-plant mass at the end of the time step
	
	Y_thetat(theta)		initial whole-plant mass
	
	s_t					starting seeding level
	
	s_0					starting seeding level used to make sure we don't have issues of derivatives getting mixed up by double-counting
	
	s_Deltat			ending seeding level
	
	diurnal_t			stores diurnal status at time t
	
	diurnal(n)			stores diurnal status at each runge-kutta step 
	
	diurnal_Deltat		stores diurnal status at time t plus delta t
	
	past_0_7			stores if past the stage 0.7 used to determine when to turn on photosynthesis because leaves are now open and free of cotyleadon
	
	past_1_2			stores if past the stage 1.2 used to determine when to turn off uptake from seed stores 
	
	/*parameters which store the final conditions*/
	x_leaf				fraction of biomass that is the leaf dependent on seeding level

	x_root				fraction of biomass that is the leaf dependent on seeding level

	x_seed				fraction of biomass that is the leaf dependent on seeding level
	
	x_stem				fraction of biomass that is the leaf dependent on seeding level
	
	s_tm1				seeding level on a scale of 0 (no seeds) to 1 max seeds wt
	
	/*parameters storing previous parameter values to allow calculations of derivatives particularly across sharp boundaries like day and night */
	t_m2				time before t_m1
	
	t_m1				time before previous hours
	
	s_tm1				stores previous seeding level
	
	s_tm2				stores seeding level prior to the previous level

	leafgrowth_ld		stores the leaf growth at the last recorded dark point	
	
	leafgrowth_ldm1		stores to leaf growth at the second to last recorded dark point
	
	leafgrowth_dl		stores the leaf growth at the last recorded light point
	
	leafgrowth_dlm1		stores the leaf growth at the second to last recorded light point
	
	omega_tm1			stores previous value of omega
	
	omega_tm2			stores previous to previous value of omega
	
	eta_tm1				stores previous value of eta
	
	eta_tm2				stores previous to previous value of eta
	
	lambda_tm1			stores previous value of lambda
	
	lambda_tm2			stores previous to previous value of lambda
	
	/*parameters updated through Runge Kutta steps*/
	t_n					current timestep
	
	Y_n(N)				keeps track to mass for the partial steps
	
	Y_thetan(N,theta)	keeps track of tissue mass for the partial steps
	
	k(N)				stores derivative estimate calculated at each Runge Kutta step

	omega_t				stores how the root tissue mass fraction changes
	
	eta_t				stores how the seed tissue mass fraction changes
	
	lambda_t			stores how the stem tissue mass fraction changes

	zeta_t				stores how the root tissue mass fraction changes in a more complex way
	
	rho_t				stores how the seed tissue mass fraction changes in a more complex way
	
	delta_t				stores how the stem tissue mass fraction changes in a more complex way
	
	Gamma_t(j)			performs multiple application Trapazoidal rule integration on metabolite concentration

	v_n(N,J)			stores reaction rate history for this timestep for trapezoidal integration

	curr_step(n)		value of 1 if the current Runge-Kutta step value of 0 otherwise

	xi_t				stores a parameter used to make formulations simpler

	/*derivatives which need to get calculated to determine k(n)*/
	ds_dt				deruvatuve if seeding level wrt time
	
	dmu_l_dt			derivative of leaf growth rate wrt time

	dlnomegat_dt		derivative of omega wrt time
	
	dlnetat_dt			derivative of eta wrt time
	
	dlnlambdat_dt		derivative of lambda wrt time
	
	lstarchrate			rate of starch storage or uptake in the leaf
	
	ststarchrate		rate of starch storage or uptake in the stem
	
	stsucroserate		rate of sucrose storate or uptake in the stem
	
	sunrise_day			sunrise for the day in question
	
	sunset_day			sunset for the day in question
	
	/*metabolite concentration parameters*/
	Z(i,day,hour)			concentration of metabolite i at indicated day and time

	Gamma(i,day,hour)		used to calculate the the trapezoid rule-based integral
	
	/*PARAMETER USED FOR FVA*/
	sigma(j)								parameter to select individual reaction for FVA
	
	/*stores when to terminate loops because maximum time is reached*/
	stop
	
;

VARIABLES

	v(j)		reaction rate
	obj			objective variables
	
;

*****************************************************************************
****************************** Symbols Defined ******************************
*****************************************************************************

*define symbols from the specifications file
c_theta('leaf') = specs('c_leaf');
c_theta('root') = specs('c_root');
c_theta('seed') = specs('c_seed');
c_theta('stem') = specs('c_stem');

x_leafinit = specs('x_leafinit');
x_rootinit = specs('x_rootinit');
x_seedinit = specs('x_seedinit');
x_steminit = specs('x_steminit');

stage_0_7 = specs('stage_0_7');
stage_1_2 = specs('stage_1_2');
Flowering_Start = specs('Flowering_Start');
Flowering_End = specs('Flowering_End');
Seed_Loss_Start = specs('Seed_Loss_Start');
Seed_Loss_End = specs('Seed_Loss_End');
Max_Time = specs('Max_Time');

inf_const = specs('inf_const');

light_up = specs('light_up');
A_l = specs('A_l');
f_l = specs('f_l');
b_l = specs('b_l');
A_st1 = specs('A_st1');
f_st1 = specs('f_st1');
b_st1 = specs('b_st1');
A_st2 = specs('A_st2');
f_st2 = specs('f_st2');
b_st2 = specs('b_st2');
k_m = specs('k_m');
k_s = specs('k_s');
Y_G = specs('Y_G');
k_m_dying = specs('k_m_dying');
k_s_dying = specs('k_s_dying');
Seeding_Gain = specs('Seeding_Gain');
Seeding_Loss = specs('Seeding_Loss');
Y_0 = specs('Y_0');
M = specs('M');
epsilon = specs('epsilon');
tau_leaf = specs('tau_leaf');
CO2_upmax = specs('CO2_upmax');
O2_upmax = specs('O2_upmax');

*define the initial state parameters
x_t('leaf') = x_leafinit;
x_t('root') = x_rootinit;
x_t('seed') = x_seedinit;
x_t('stem') = x_steminit;

Y_thetat('leaf') = Y_0 * x_t('leaf');
Y_thetat('root') = Y_0 * x_t('root');
Y_thetat('seed') = Y_0 * x_t('seed');
Y_thetat('stem') = Y_0 * x_t('stem');

t_0 = 0;
t_Deltat = t_0 + 1;
t_m1 = t_0 - 1;
t_m2 = t_m1 - 1;

s_Deltat = (0)$((t_Deltat < Flowering_Start) OR (t_Deltat ge Seed_Loss_End)) + (Seeding_Gain * (t_Deltat - Flowering_Start + 1))$((t_Deltat ge Flowering_Start) AND (t_Deltat < Seed_Loss_Start)) + (Seeding_Gain * (t_Deltat - Flowering_Start + 1) - Seeding_Loss * (t_Deltat - Seed_Loss_Start + 1))$((t_Deltat ge Flowering_Start) AND (t_Deltat ge Seed_Loss_Start) AND (t_Deltat < Flowering_End)) + (1 - Seeding_Loss * (t_Deltat - Seed_Loss_Start + 1))$((t_Deltat ge Flowering_End) AND (t_Deltat ge Seed_Loss_Start) AND (t_Deltat < Seed_Loss_End)) + 0;
s_t = (0)$((t_0 < Flowering_Start) OR (t_0 ge Seed_Loss_End)) + (Seeding_Gain * (t_0 - Flowering_Start + 1))$((t_0 ge Flowering_Start) AND (t_0 < Seed_Loss_Start)) + (Seeding_Gain * (t_0 - Flowering_Start + 1) - Seeding_Loss * (t_0 - Seed_Loss_Start + 1))$((t_0 ge Flowering_Start) AND (t_0 ge Seed_Loss_Start) AND (t_0 < Flowering_End)) + (1 - Seeding_Loss * (t_0 - Seed_Loss_Start + 1))$((t_0 ge Flowering_End) AND (t_0 ge Seed_Loss_Start) AND (t_0 < Seed_Loss_End)) + 0;
s_tm1 = (0)$((t_m1 < Flowering_Start) OR (t_m1 ge Seed_Loss_End)) + (Seeding_Gain * (t_m1 - Flowering_Start + 1))$((t_m1 ge Flowering_Start) AND (t_m1 < Seed_Loss_Start)) + (Seeding_Gain * (t_m1 - Flowering_Start + 1) - Seeding_Loss * (t_m1 - Seed_Loss_Start + 1))$((t_m1 ge Flowering_Start) AND (t_m1 ge Seed_Loss_Start) AND (t_m1 < Flowering_End)) + (1 - Seeding_Loss * (t_m1 - Seed_Loss_Start + 1))$((t_m1 ge Flowering_End) AND (t_m1 ge Seed_Loss_Start) AND (t_m1 < Seed_Loss_End)) + 0;
s_tm2 = (0)$((t_m2 < Flowering_Start) OR (t_m2 ge Seed_Loss_End)) + (Seeding_Gain * (t_m2 - Flowering_Start + 1))$((t_m2 ge Flowering_Start) AND (t_m2 < Seed_Loss_Start)) + (Seeding_Gain * (t_m2 - Flowering_Start + 1) - Seeding_Loss * (t_m2 - Seed_Loss_Start + 1))$((t_m2 ge Flowering_Start) AND (t_m2 ge Seed_Loss_Start) AND (t_m2 < Flowering_End)) + (1 - Seeding_Loss * (t_m2 - Seed_Loss_Start + 1))$((t_m2 ge Flowering_End) AND (t_m2 ge Seed_Loss_Start) AND (t_m2 < Seed_Loss_End)) + 0;

omega_tm2 = 1;
omega_tm1 = 1;
omega_t = (x_t('root') * Y_thetat('leaf')) / (x_t('leaf') * Y_thetat('root'));

eta_tm2 = 0;
eta_tm1 = 0;
eta_t = ((x_t('seed') * Y_thetat('leaf')) / (x_t('leaf') * Y_thetat('seed')))$(Y_thetat('seed') ne 0) + (0)$(Y_thetat('seed') eq 0);

lambda_tm2 = 1;
lambda_tm1 = 1;
lambda_t = (x_t('stem') * Y_thetat('leaf')) / (x_t('leaf') * Y_thetat('stem'));

*right now assumed that for the last several hours had no growth
dmu_l_dt = 0;

*set initial metabolite concentration
Z(I,day,hour) = 0;
Z_0(I) = 0;

*set bounds of the reactions
LB(J)$(rxntype(J) eq -1) = -M; 
UB(J)$(rxntype(J) eq -1) = 0; 

LB(J)$(rxntype(J) eq 0) = -M; 
UB(J)$(rxntype(J) eq 0) = M; 

LB(J)$(rxntype(J) eq 1) = 0; 
UB(J)$(rxntype(J) eq 1) = M;

LB(reg) = 0;
UB(reg) = 0;

*set physiological limits for carbon dioxide and water uptake
LB('R164ex[l]') = - CO2_upmax;
LB('R161ex[l]') = - O2_upmax;

LB('R044ex[r]') = - 0.1 * O2_upmax;
LB('R116ex[se]') = - 0.1 * O2_upmax;
LB('R059ex[st]') = - 0.1 * O2_upmax;

*allow loops to run
stop = 0;

*****************************************************************************
*************************** Equations Initialized ***************************
*****************************************************************************

EQUATIONS

	/* Outer objective problem*/
	growthobj         			Objective function function
	nleffobj					linear photonic efficiency and seed growth objective
	growthpFAobj				leaf growth objective seed fatty acid storage objective
	lineffpFAobj				leaf photonic efficiency and seed fatty acid storage objective
	
	/* leaf tissue primal and dual*/
	mass_balance(i)				Mass balance for metabolites allowing storage
	mass_balance_2(i)			Limit rate of storage at each point
	mass_balance_3(i)			strict mass balance for certain metabolites
	
	/*constraints to enforce reaction bounds*/
	ub_const(J)
	lb_const(J)
	
	/*tissue interaction constraints. Naming format:*/
	/*Compound _ tissue export uptake _ tissue export uptake _ const*/
	/*the constraint relating the  export (e) of water in the root to the uptake of water in the root (e)*/
	/*multiple tissues are listed in alphabetical order where the underscore between tissues is where the equality/inequality is*/
	/*water constraints*/
	C00001_stem_root_const
	C00001_leaf_seed_stem_const
	
	/*oxygen constraints*/
	net_C00007_const
	C00007_leaf_leaf1
	C00007_leaf_leaf2
	C00007_root
	C00007_seed
	C00007_stem
	
	/*phosphate constraints*/
	C00009_stem_root_const
	C00009_stem_stem_const
	C00009_leaf_seed_stem_const
	
	/*second the plant must net uptake carbon dioxide*/
	net_C00011_const
	C00011_leaf_leaf1
	C00011_leaf_leaf2
	
	/*nitrate constraints*/ 
	C00244_root_root_const
	C00244_stem_root_const
	C00244_leaf_seed_stem_const
	
	/*L-glutamate constraints*/
	C00025_leaf_stem_const
	C00025_stem_stem_const
	C00025_seed_stem_const
	
	/*L-glutamate constraints*/
	C00037_leaf_stem_const
	C00037_stem_stem_const
	C00037_seed_stem_const
	
	/*L-alanine constraints*/
	C00041_leaf_stem_const
	C00041_stem_stem_const
	C00041_seed_stem_const
	
	/*L-lysine constraints*/
	C00047_leaf_stem_const
	C00047_stem_stem_const
	C00047_seed_stem_const
	
	/*L-Aspartate constraints*/
	C00049_leaf_seed_stem_const
	C00049_stem_stem_const
	C00049_root_stem_const
	
	/*sulfate constraints*/
	C00059_stem_root_const
	C00059_stem_stem_const
	C00059_leaf_seed_stem_const
	
	/*l-arginine constraints*/ 
	C00062_leaf_stem_const
	C00062_stem_stem_const
	C00062_seed_stem_const
	
	/*l-glutamine constraints*/
	C00064_leaf_seed_stem_const
	C00064_stem_stem_const
	C00064_stem_root_const
	
	/*l-serine constraints*/
	C00065_leaf_stem_const
	C00065_stem_stem_const
	C00065_seed_stem_const
	
	/*l-methionine constraints*/
	C00073_leaf_stem_const
	C00073_stem_stem_const
	C00073_seed_stem_const
	
	/*l-tryptophan constraints*/
	C00078_stem_leaf_const
	C00078_stem_stem_const
	C00078_seed_stem_const
	
	/*proton constraints*/
	C00080_root_stem_const
	C00080_stem_seed_leaf_const
	
	/*sucrose constraints*/
	C00089_stem_leaf_const
	C00089_stem_stem_const
	C00089_root_seed_stem_const
	
	/*L-cysteine constraints*/
	C00097_leaf_stem_const
	C00097_stem_stem_const
	C00097_seed_stem_const
	
	/*L-proline constraints*/
	C00148_leaf_stem_const
	C00148_stem_stem_const
	C00148_seed_stem_const
	
	/*L-valine constraints*/
	C00183_leaf_stem_const
	C00183_stem_stem_const
	C00183_seed_stem_const
	
	/*L-threonine constraints*/
	C00188_leaf_stem_const
	C00188_stem_stem_const
	C00188_seed_stem_const
	
	/*light constraints*/
	C00205_leaf
	
	/*starch constraints */
	C00369_leaf_1
	C00369_leaf_2
	C00369_leaf_3
	C00369_leaf_4
	C00369_leaf_5
	C00369_stem
	C00089_stem
	
	/*uptake of seed stores for the growth of the seedling*/
	up_seed_const(seed_stores)
	
	/*biomass ratios constraining growth*/
	bio_leaf_root
	bio_leaf_stem
	bio_leaf_seed
	
	/*constraints to force respiration*/
	transpiration_const
	
	/*constraint to force biomass_loss = biomass_wt when seed is loosing mass*/
	seed_mass_loss
	
	/*equation for performing FVA on the the model based on the growth rate*/
	FVA_obj
	
	/*scenescenece and maintenance costs*/
	s_and_m_cost_leaf
	s_and_m_cost_root
	s_and_m_cost_seed
	s_and_m_cost_stem

;

*****************************************************************************
***************************** Equations Defined *****************************
*****************************************************************************

*Outer objective problem
growthobj..			obj =e= v('lbiomasssi') + v('rbiomasssi') + v('stbiomasssi') - epsilon * epsilon * sum(j,v(j));

*this formulation should ensure some ability to make use of stored metaboltes
mass_balance(I)..	sum(J, S(I,J) * v(J)) =g= -Z_0(I);
mass_balance_2(I)..	sum(J, S(I,J) * v(J)) =l= 10;
mass_balance_3(I_prime)..	sum(J, S(I_prime,J) * v(J)) =e= 0;
*mass_balance_3(I)..	sum(J, S(I,J) * v(J)) =e= 0;

*reaction bounds constraint
ub_const(J)..		v(J) =l= UB(J);
lb_const(J)..		v(J) =g= LB(J);

*community interaction constraints
*water (C00001)
*In root R040ex accounts for uptake of water from soil
C00001_stem_root_const..		(-v('R063ex[st]')) * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf')) =e= v('R039ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root'));							!!water uptake by stem must equal to that exported by root since the stem is the only destination
C00001_leaf_seed_stem_const..	(-v('R115ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) + (-v('R160ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'))) =e= v('R062ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));	!!water exported from stem has 2 destinations possible: leaf or seed both of which uptake water

*oxygen (C00007) 
net_C00007_const..			(v('R044ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root')) + v('R116ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed')) + v('R059ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) + v('R161ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'))) * (1 - 2 * diurnal_t) =l= 0;		!!amount of oxygen output by plant in light must be non-negative
C00007_leaf_leaf1..			v('R161ex[l]') =g= -O2_upmax * (1/sum(N1,curr_step(N1) * x_n(N1,'leaf'))) * (1 - diurnal_t);
C00007_leaf_leaf2..			v('R161ex[l]') =l= M * diurnal_t;
*these just scale the bounds of oxygen uptake of other tissues by an order of magnitude
*since assuming leaf much better than other tissues
*also no diurnal switch on gas direction in these tissues
C00007_root..			v('R044ex[r]') =g= -0.1 * O2_upmax * (1/sum(N1,curr_step(N1) * x_n(N1,'leaf')));
C00007_seed..			v('R116ex[se]') =g= -0.1 * O2_upmax * (1/sum(N1,curr_step(N1) * x_n(N1,'leaf')));
C00007_stem..			v('R059ex[st]') =g= -0.1 * O2_upmax * (1/sum(N1,curr_step(N1) * x_n(N1,'leaf')));

*phosphate (C00009)
C00009_stem_root_const..		(-v('R065ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R041ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root'));							!!export of phosphate from root goes to the stem
C00009_stem_stem_const..		v('R064ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =l= (-v('R065ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));							!!stem exports less phosphate than it uptakes (stem uses phosphate)
C00009_leaf_seed_stem_const..	(-v('R162ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'))) + (-v('R117ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R064ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));	!!phosphate exported from stem either goes to seed or leaf

*carbon dioxide (C00011) 
net_C00011_const..			(v('R164ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf')) + v('R066ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) + v('R125ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed')) + v('R043ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root'))) * (1 - 2 * diurnal_t) =g= 0;					!!uptake carbon dioxide must be non-negative
C00011_leaf_leaf1..			v('R164ex[l]') =g= -CO2_upmax  * (1/sum(N1,curr_step(N1) * x_n(N1,'leaf'))) * diurnal_t;
C00011_leaf_leaf2..			v('R164ex[l]') =l= M * (1 - diurnal_t);

*ammonia (C00244)
C00244_root_root_const..		v('R038ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root')) =l= (-v('R037ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root')));							!!root uses ammonia exports less than it takes up
C00244_stem_root_const..		(-v('R060ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R038ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root'));						!!all ammonia exported by root goes to stem
C00244_leaf_seed_stem_const..	(-v('R157ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'))) + (-v('R113ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R061ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));	!!ammonia exported by stems is either taken up the the leaf or by the seed

*l-glutamate (amino acid) (C00025)
C00025_leaf_stem_const..		(-v('R075ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R177ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));						!!l-alanine export from leaf where it is synthesized to the stem for transport
C00025_stem_stem_const..		(-v('R075ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =g= v('R076ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));						!!l-alanine export from leaf where it is synthesized to the stem for transport
C00025_seed_stem_const..		(-v('R130ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R076ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-alanine is imported into the seed for use

*l-glycine (amino acid) (C00037)
C00037_leaf_stem_const..		(-v('R077ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R178ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));						!!l-alanine export from leaf where it is synthesized to the stem for transport
C00037_stem_stem_const..		(-v('R077ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R078ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-alanine export from leaf where it is synthesized to the stem for transport
C00037_seed_stem_const..		(-v('R131ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R078ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-alanine is imported into the seed for use

*L-alanine (amino acid) (C00041)
C00041_leaf_stem_const..		(-v('R053ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R173ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));						!!l-alanine export from leaf where it is synthesized to the stem for transport
C00041_stem_stem_const..		(-v('R053ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =g= v('R052ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-alanine export from leaf where it is synthesized to the stem for transport
C00041_seed_stem_const..		(-v('R112ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R052ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-alanine is imported into the seed for use

*L-lysine (amino acid) (C00047)
C00047_leaf_stem_const..		(-v('R045ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R169ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));						!!l-lysine export from leaf where it is synthesized to the stem for transport
C00047_stem_stem_const..		(-v('R045ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =g= v('R044ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-lysine export from leaf where it is synthesized to the stem for transport
C00047_seed_stem_const..		(-v('R121ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R044ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-lysine is imported into the seed for use

*L-aspartate (amino acid) (C00049)
C00049_leaf_seed_stem_const..	(-v('R069ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R174ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf')) + v('R127ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'));						!!l-lysine export from leaf where it is synthesized to the stem for transport
C00049_stem_stem_const..		(-v('R069ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =g= v('R070ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-lysine is imported into the seed for use
C00049_root_stem_const..		(v('R051ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root'))) =e= -v('R069ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-lysine is imported into the seed for use

*Sulfate (C00059)
C00059_stem_root_const..		(-v('R058ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R045ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root'));						!!sulfate transfered from root to stem
C00059_stem_stem_const..		v('R057ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =l= (-v('R058ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));						!!all sulfate uptaken by stem is exported by stem (stem doesn't use sulfate)
C00059_leaf_seed_stem_const..	(-v('R158ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'))) + (-v('R114ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R057ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));	!!all sulfate exported by stem either goes to the seed or the the leaf

*l-arginine (amino acid) (C00062)
C00062_leaf_stem_const..		(-v('R041ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R165ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));						!!l-arginine export from leaf where it is synthesized to the stem for transport
C00062_stem_stem_const..		v('R040ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =e= (-v('R041ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));						!!l-arginine export from stem equal to import because stem has no AA metabolism will not use any AA
C00062_seed_stem_const..		(-v('R119ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R040ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-arginine is imported into the seed for use

*l-glutamine (amino acid) (C00064)
C00064_stem_stem_const..		v('R038ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =l= (-v('R039ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));						!!l-glutamine export from stem equal to import because stem has no AA metabolism will not use any AA
C00064_leaf_seed_stem_const..	(-v('R110ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) + (-v('R156ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'))) =e= v('R038ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-glutamine is imported into the seed for use
C00064_stem_root_const..		(-v('R039ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R052ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root'));

*l-serine (amino acid) (C00065)
C00065_leaf_stem_const..		(-v('R073ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R176ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));						!!l-glutamine export from leaf where it is synthesized to the stem for transport
C00065_stem_stem_const..		v('R074ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =l= (-v('R073ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));						!!l-glutamine export from stem equal to import because stem has no AA metabolism will not use any AA
C00065_seed_stem_const..		(-v('R129ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R074ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-glutamine is imported into the seed for use

*l-methionine (amino acid) (C00073)
C00073_leaf_stem_const..		(-v('R047ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R166ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));						!!l-methionine export from leaf where it is synthesized to the stem for transport
C00073_stem_stem_const..		v('R046ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =l= (-v('R047ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));						!!l-methionine export from stem equal to import because stem has no AA metabolism will not use any AA
C00073_seed_stem_const..		(-v('R122ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R046ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-methionine is imported into the seed for use

*l-tryptophan (amino acid) (C00078)
C00078_stem_leaf_const..		(-v('R051ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R168ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));						!!l-tryptophan export from leaf where it is synthesized to the stem for transport
C00078_stem_stem_const..		v('R050ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =l= (-v('R051ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));						!!l-tryptophan export from stem equal to import because stem has no AA metabolism will not use any AA
C00078_seed_stem_const..		(-v('R124ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R050ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));						!!l-tryptophan is imported into the seed for use

*proton constraints
*note that root uptake or export from soil will be handled by R049ex[r]
*note that stem uptake or export from root will be handled by R067ex[st]
C00080_root_stem_const..		v('R050ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root')) =e= -v('R067ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));
C00080_stem_seed_leaf_const..	(v('R171ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'))) + (v('R126ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= -v('R068ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));	!!all protons exported by stem either goes to the seed or the the leaf

*Sucrose (C00089)
C00089_stem_leaf_const..		(-v('R054ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R172ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));						!!sucrose export from leaf is stem uptake
C00089_stem_stem_const..		v('R055ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =l= (-v('R054ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));						!!sucrose export from stem less than uptake from stem (stem uses sucrose)
C00089_root_seed_stem_const..	(-v('R111ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) + (-v('R047ex[r]') * sum(N1, curr_step(N1) * Y_thetan(N1,'root'))) =e= v('R055ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));	!!exported sucrose from stem goes to seed or root

*L-cysteine (C00097)
C00097_leaf_stem_const..		(-v('R043ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R170ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));		!!l-cysteine export from leaf where it is synthesized to the stem for transport
C00097_stem_stem_const..		v('R042ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =l= (-v('R043ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));		!!l-cysteine export from stem equal to import because stem has no AA metabolism will not use any AA
C00097_seed_stem_const..		(-v('R120ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R042ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));		!!l-cysteine is imported into the seed for use

*L-proline (C00148)
C00148_leaf_stem_const..		(-v('R049ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R167ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));		!!l-proline export from leaf where it is synthesized to the stem for transport
C00148_stem_stem_const..		v('R048ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =l= (-v('R049ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));		!!l-proline export from stem equal to import because stem has no AA metabolism will not use any AA
C00148_seed_stem_const..		(-v('R123ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R048ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));		!!l-proline is imported into the seed for use

*L-valine (C00183)
C00183_leaf_stem_const..		(-v('R079ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R179ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));		!!l-proline export from leaf where it is synthesized to the stem for transport
C00183_stem_stem_const..		v('R080ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =e= (-v('R079ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));		!!l-proline export from stem equal to import because stem has no AA metabolism will not use any AA
C00183_seed_stem_const..		(-v('R132ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R080ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));		!!l-proline is imported into the seed for use

*L-threonine (C00188)
C00188_leaf_stem_const..		(-v('R071ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'))) =e= v('R175ex[l]') * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));		!!l-proline export from leaf where it is synthesized to the stem for transport
C00188_stem_stem_const..		v('R072ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')) =l= (-v('R071ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem')));		!!l-proline export from stem equal to import because stem has no AA metabolism will not use any AA
C00188_seed_stem_const..		(-v('R128ex[se]') * sum(N1, curr_step(N1) * Y_thetan(N1,'seed'))) =e= v('R072ex[st]') * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));		!!l-proline is imported into the seed for use

*light constraints
C00205_leaf..					v('R163ex[l]') =g= (-light_up * (1/sum(N1,curr_step(N1) * x_n(N1,'leaf')))) * diurnal_t;

*constraints on storage and uptake of carbohydrates in leaf and stem
*no diurnal patterns of carbohydrate storage if the seedling cannot use light for photosynthesis
C00369_leaf_1..					v('R151si[l]') + v('R152si[l]') =e= (A_l * sin(f_l * (t_n + b_l))) * past_0_7;
C00369_leaf_2..					v('R152si[l]') =g= -M * past_0_7 + (M * past_0_7)$((A_l * sin(f_l * (t_n + b_l))) > 0);
C00369_leaf_3..					v('R151si[l]') =g= -M * past_0_7 + (M * past_0_7)$((A_l * sin(f_l * (t_n + b_l))) > 0);
C00369_leaf_4..					v('R152si[l]') =l= M * past_0_7 - (M * past_0_7)$((A_l * sin(f_l * (t_n + b_l))) < 0);
C00369_leaf_5..					v('R151si[l]') =l= M * past_0_7 - (M * past_0_7)$((A_l * sin(f_l * (t_n + b_l))) < 0);
C00369_stem..					v('R034si[st]') =e= (A_st1 * sin(f_st1 * (t_n + b_st1))) * past_0_7;
C00089_stem..					v('R037si[st]') =e= (A_st2 * sin(f_st2 * (t_n + b_st2))) * past_0_7;

*constraints on the uptake of metabolites from the seed stores to feed the seedling
up_seed_const(seed_stores)..	v(seed_stores) =g= (- (1 / sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'))) * up_rate(seed_stores)) * (1 - past_1_2) + 0;

*biomass ratios constraining growth
bio_leaf_root..					v('rbiomasssi') =e= log(omega_t) + v('lbiomasssi');
bio_leaf_stem..					v('sebiomasssi') =e= (0)$((sum(N1, curr_step(N1) * Y_thetan(N1,'seed')) eq 0) and (s_t eq 0)) + (log(eta_t) + v('lbiomasssi'))$((eta_t ne 0) AND (s_t ne 0) AND (sum(N1, curr_step(N1) * Y_thetan(N1,'seed')) ne 0)) + (-v('lbiomasssi'))$((sum(N1, curr_step(N1) * Y_thetan(N1,'seed')) ne 0) and (s_t eq 0)) + (v('lbiomasssi'))$((sum(N1, curr_step(N1) * Y_thetan(N1,'seed')) eq 0) and (s_t ne 0));
bio_leaf_seed..					v('stbiomasssi') =e= log(lambda_t) + v('lbiomasssi');

*constraints to force respiration
transpiration_const..			v('R159ex[l]') =e= tau_leaf  * (1/sum(N1,curr_step(N1) * x_n(N1,'leaf'))) * diurnal_t;

*constraint to force biomass_loss = biomass_wt when seed is loosing mass
seed_mass_loss..				v('sebio_loss_ws') =e= (- v('sebiomass_ws'))$((t_n ge Flowering_End) AND (t_n ge Seed_Loss_Start)) + (0)$(t_n < Seed_Loss_Start);
	
*equation for performing FVA on the the model based on the growth rate
FVA_obj..						obj =e= sum(j, sigma(j) * v(j));
	
*scenescenece and maintenance costs
s_and_m_cost_leaf..				v('lbiomass_scen') =e= ((k_m + k_s)/24) * sum(N1, curr_step(N1) * Y_thetan(N1,'leaf'));
s_and_m_cost_root..				v('rbiomass_scen') =e= ((k_m + k_s)/24) * sum(N1, curr_step(N1) * Y_thetan(N1,'root'));
s_and_m_cost_seed..				v('sebiomass_scen') =e= (((k_m + k_s)/24) * sum(N1, curr_step(N1) * Y_thetan(N1,'seed')))$((t_n > Flowering_Start) AND (t_n < Flowering_End));
s_and_m_cost_stem..				v('stbiomass_scen') =e= ((k_m + k_s)/24) * sum(N1, curr_step(N1) * Y_thetan(N1,'stem'));

*****************************************************************************
******************************* Models Defined ******************************
*****************************************************************************
MODEL PATH780
/

	growthobj
	mass_balance
	mass_balance_2
	mass_balance_3
	ub_const
	lb_const
	C00001_stem_root_const
	C00001_leaf_seed_stem_const
	net_C00007_const
	C00007_leaf_leaf1
	C00007_leaf_leaf2
	C00007_root
	C00007_seed
	C00007_stem
	C00009_stem_root_const
	C00009_stem_stem_const
	C00009_leaf_seed_stem_const
	net_C00011_const
	C00011_leaf_leaf1
	C00011_leaf_leaf2
	C00244_root_root_const
	C00244_stem_root_const
	C00244_leaf_seed_stem_const
	C00025_leaf_stem_const
	C00025_stem_stem_const
	C00025_seed_stem_const
	C00037_leaf_stem_const
	C00037_stem_stem_const
	C00037_seed_stem_const
	C00041_leaf_stem_const
	C00041_stem_stem_const
	C00041_seed_stem_const
	C00047_leaf_stem_const
	C00047_stem_stem_const
	C00047_seed_stem_const
	C00049_leaf_seed_stem_const
	C00049_stem_stem_const
	C00049_root_stem_const
	C00059_stem_root_const
	C00059_stem_stem_const
	C00059_leaf_seed_stem_const
	C00062_leaf_stem_const
	C00062_stem_stem_const
	C00062_seed_stem_const
	C00064_leaf_seed_stem_const
	C00064_stem_stem_const
	C00064_stem_root_const
	C00065_leaf_stem_const
	C00065_stem_stem_const
	C00065_seed_stem_const
	C00073_leaf_stem_const
	C00073_stem_stem_const
	C00073_seed_stem_const
	C00078_stem_leaf_const
	C00078_stem_stem_const
	C00078_seed_stem_const
	C00080_root_stem_const
	C00080_stem_seed_leaf_const
	C00089_stem_leaf_const
	C00089_stem_stem_const
	C00089_root_seed_stem_const
	C00097_leaf_stem_const
	C00097_stem_stem_const
	C00097_seed_stem_const
	C00148_leaf_stem_const
	C00148_stem_stem_const
	C00148_seed_stem_const
	C00183_leaf_stem_const
	C00183_stem_stem_const
	C00183_seed_stem_const
	C00188_leaf_stem_const
	C00188_stem_stem_const
	C00188_seed_stem_const
	C00205_leaf
	C00369_leaf_1
	C00369_leaf_2
	C00369_leaf_3
	C00369_leaf_4
	C00369_leaf_5
	C00369_stem
	C00089_stem
	up_seed_const
	bio_leaf_root
	bio_leaf_seed
	bio_leaf_stem
	transpiration_const
	seed_mass_loss
	s_and_m_cost_leaf
	s_and_m_cost_root
	s_and_m_cost_seed
	s_and_m_cost_stem

/
;

MODEL FVA_PATH780
/

	FVA_obj
	mass_balance
	mass_balance_2
	mass_balance_3
	ub_const
	lb_const
	C00001_stem_root_const
	C00001_leaf_seed_stem_const
	net_C00007_const
	C00007_leaf_leaf1
	C00007_leaf_leaf2
	C00007_root
	C00007_seed
	C00007_stem
	C00009_stem_root_const
	C00009_stem_stem_const
	C00009_leaf_seed_stem_const
	net_C00011_const
	C00011_leaf_leaf1
	C00011_leaf_leaf2
	C00244_root_root_const
	C00244_stem_root_const
	C00244_leaf_seed_stem_const
	C00025_leaf_stem_const
	C00025_stem_stem_const
	C00025_seed_stem_const
	C00037_leaf_stem_const
	C00037_stem_stem_const
	C00037_seed_stem_const
	C00041_leaf_stem_const
	C00041_stem_stem_const
	C00041_seed_stem_const
	C00047_leaf_stem_const
	C00047_stem_stem_const
	C00047_seed_stem_const
	C00049_leaf_seed_stem_const
	C00049_stem_stem_const
	C00049_root_stem_const
	C00059_stem_root_const
	C00059_stem_stem_const
	C00059_leaf_seed_stem_const
	C00062_leaf_stem_const
	C00062_stem_stem_const
	C00062_seed_stem_const
	C00064_leaf_seed_stem_const
	C00064_stem_stem_const
	C00064_stem_root_const
	C00065_leaf_stem_const
	C00065_stem_stem_const
	C00065_seed_stem_const
	C00073_leaf_stem_const
	C00073_stem_stem_const
	C00073_seed_stem_const
	C00078_stem_leaf_const
	C00078_stem_stem_const
	C00078_seed_stem_const
	C00080_root_stem_const
	C00080_stem_seed_leaf_const
	C00089_stem_leaf_const
	C00089_stem_stem_const
	C00089_root_seed_stem_const
	C00097_leaf_stem_const
	C00097_stem_stem_const
	C00097_seed_stem_const
	C00148_leaf_stem_const
	C00148_stem_stem_const
	C00148_seed_stem_const
	C00183_leaf_stem_const
	C00183_stem_stem_const
	C00183_seed_stem_const
	C00188_leaf_stem_const
	C00188_stem_stem_const
	C00188_seed_stem_const
	C00205_leaf
	C00369_leaf_1
	C00369_leaf_2
	C00369_leaf_3
	C00369_leaf_4
	C00369_leaf_5
	C00369_stem
	C00089_stem
	up_seed_const
	bio_leaf_root
	bio_leaf_seed
	bio_leaf_stem
	transpiration_const
	seed_mass_loss
	s_and_m_cost_leaf
	s_and_m_cost_root
	s_and_m_cost_seed
	s_and_m_cost_stem

/
;

PATH780.optfile = 1;
FVA_PATH780.optfile = 1;

*****************************************************************************
**************************** PREPARE OUTPUT FILES ***************************
*****************************************************************************
*print FBA-based optimal flux results
*prints the rate of each reaction at a specified date and hour
FILE FLUXSHEET /Flux_sheet.csv/;
PUT FLUXSHEET;
FLUXSHEET.pc=5;
FLUXSHEET.pw=32767;
FLUXSHEET.nr=2;

PUT "day","hour","t_0";

LOOP(j1,

	PUT j1.tl;

);

PUT /;

PUTCLOSE;

*print CBA-based optimal flux results
*prints the rate of each reaction at a specified date and hour
FILE CONCSHEET /Concentration_sheet.csv/;
PUT CONCSHEET;
CONCSHEET.pc=5;
CONCSHEET.pw=32767;
CONCSHEET.nr=2;

PUT "day","hour","t_0";

LOOP(i1,

	PUT i1.tl;

);

PUT /;

PUTCLOSE;

*create a file for troubleshooting
FILE TROUBLESHOOT /troubleshooting_revised.txt/;
PUT TROUBLESHOOT;
TROUBLESHOOT.nr = 2;
PUT "Troubleshooting file"/;
PUT "Started Date: ",system.date,system.tab,"Time: ", system.time;

PUTCLOSE;

FILE MASSSHEET /mass_sheet.csv/;
PUT MASSSHEET;
PUT "Plant and tissue mass output file"/;
MASSSHEET.pc=5;
PUT "Day","Hour","t_0","Y_0","Y_leaf0","Y_root0","Y_seed0","Y_stem0,";
MASSSHEET.pc=1;

LOOP(N,

	PUT "Y_n(",n.tl,"),";

);

MASSSHEET.pc=5;
PUT "dY_dt|rk est"/;

PUTCLOSE;

*****************************************************************************
****************************** SOLVE THE MODELS *****************************
*****************************************************************************

LOOP(day$(stop eq 0),

	sunrise_day = sunrise(day);
	sunset_day = sunset(day);

	LOOP(hour$(stop eq 0), 
		
		t_Deltat = t_0 + 1;
		Z(I,day,hour) = Z_0(I);
		
		/*restart the N minus the last element set*/
		Nm1(N) = NO;
		
		/*initialize values of k and v_n apparently necessary to prevent compilation errors*/
		k(N) = 0; 
		v_n(N,J) = 0;
		Y_n(N) = 0;
		
		TROUBLESHOOT.ap = 1;
		PUT TROUBLESHOOT;
		TROUBLESHOOT.nr = 1;
		
		PUT "#####################################################################################"/;
		PUT "#################### BEGIN TIMESTEP",t_0," to ",t_Deltat," ####################"/;
		PUT "#####################################################################################"/;
		PUT /;
		
		PUTCLOSE;
		
		/*update and save tissue mass and mass fraction*/
		s_0 = (0)$((t_0 <= Flowering_Start) OR (t_0 ge Seed_Loss_End)) + (Seeding_Gain * (t_0 - Flowering_Start + 1))$((t_0 ge Flowering_Start) AND (t_0 < Seed_Loss_Start)) + (Seeding_Gain * (t_0 - Flowering_Start + 1) - Seeding_Loss * (t_0 - Seed_Loss_Start + 1))$((t_0 ge Flowering_Start) AND (t_0 ge Seed_Loss_Start) AND (t_0 < Flowering_End)) + (1 - Seeding_Loss * (t_0 - Seed_Loss_Start + 1))$((t_0 ge Flowering_End) AND (t_0 ge Seed_Loss_Start) AND (t_0 < Seed_Loss_End)) + 0;
		
		IF ((s_0 < 0),
			
			s_0 = 0;
			
		);
		
		x_t('leaf') = c_theta('leaf') * s_0 + x_leafinit;
		x_t('root') = c_theta('root') * s_0 + x_rootinit;
		x_t('seed') = c_theta('seed') * s_0 + x_seedinit;
		x_t('stem') = c_theta('stem') * s_0 + x_steminit;

		Y_thetat('leaf') = Y_0 * x_t('leaf');
		Y_thetat('root') = Y_0 * x_t('root');
		Y_thetat('seed') = Y_0 * x_t('seed');
		Y_thetat('stem') = Y_0 * x_t('stem');
		
		LOOP(N,
		
			/*Loops for the Runge Kutta steps*/
			/*update what parameters need updating for the given time*/
			curr_step(n1) = 0;
			curr_step(n) = 1;
			
			/*update the time and plant mass*/
			t_n = t_0 + c(n) * 1;
			Y_n(N) = Y_0 + 1 * (sum(Nm1, a(N,Nm1) * k(Nm1)));
			
			/*update k_m and k_s if needed*/
			IF ((t_n le Seed_Loss_Start),
			
				/*do nothing*/
			
			ELSEIF (t_n > Seed_Loss_Start),
				
				/*update maintenance and senescence*/
				k_m = k_m_dying;
				k_s = k_s_dying;
			
			);
			
			/*determine diurnal status*/
			IF(((timeOfDay(hour) >= sunrise_day) AND (timeOfDay(hour) < sunset_day) AND (t_n > stage_0_7)),
			
				/*if there is light available diurnal status will be saved as one*/
				diurnal_t = 1;
				
			ELSE
			
				/*if there is no light available, diurnal status will be saved as zero*/
				diurnal_t = 0;
			
			);
			
			IF ((t_n > stage_0_7),
			
				past_0_7 = 1;
			
			ELSE
			
				past_0_7 = 0;
			
			);
			
			IF ((t_n > stage_1_2),
			
				past_1_2 = 1;
			
			ELSE
			
				past_1_2 = 0;
			
			);
			
			/*update parameters for solving the model*/
			s_tm2 = s_tm1;
			s_tm1 = s_t;
			s_t = (0)$((t_n < Flowering_Start) OR (t_n ge Seed_Loss_End)) + (Seeding_Gain * (t_n - Flowering_Start + 1))$((t_n ge Flowering_Start) AND (t_n < Seed_Loss_Start)) + (Seeding_Gain * (t_n - Flowering_Start + 1) - Seeding_Loss * (t_n - Seed_Loss_Start + 1))$((t_n ge Flowering_Start) AND (t_n ge Seed_Loss_Start) AND (t_n < Flowering_End)) + (1 - Seeding_Loss * (t_n - Seed_Loss_Start + 1))$((t_n ge Flowering_End) AND (t_n ge Seed_Loss_Start) AND (t_n < Seed_Loss_End)) + 0;
			
			/*correction for if the slopes of loss and gain don't quite match force seeding to zero*/
			/*if seeding drops below zero*/
			IF((s_t < 0),
			
				s_t = 0;
			
			);
			
			x_n(N,'leaf') = c_theta('leaf') * s_t + x_leafinit;
			x_n(N,'root') = c_theta('root') * s_t + x_rootinit;
			x_n(N,'seed') = c_theta('seed') * s_t + x_seedinit;
			x_n(N,'stem') = c_theta('stem') * s_t + x_steminit;
			
			Y_thetan(N,'leaf') = Y_n(N) * x_n(N,'leaf');
			Y_thetan(N,'root') = Y_n(N) * x_n(N,'root');
			Y_thetan(N,'seed') = Y_n(N) * x_n(N,'seed');
			Y_thetan(N,'stem') = Y_n(N) * x_n(N,'stem');
			
			/*if no seeds turn off all seed reactions*/
			IF(((Y_thetan(N,'seed') eq 0) AND (s_t eq 0)),
			
				UB(jseed) = 0;
				LB(jseed) = 0;
			
			/*if seeds at start but not at end of time window*/
			/*turn off all reactions except loss of seed mass*/
			ELSEIF ((Y_thetan(N,'seed') ne 0) AND (s_t eq 0)),
			
				UB(jseed) = 0;
				LB(jseed) = 0;
				
				/*allow biomass to go negative*/
				LB('sebiomass_ws') = -M;
				UB('sebiomass_ws') = 0;

				/*allow positive loss of biomass*/
				LB('sebio_loss_ws') = 0;
				UB('sebio_loss_ws') = M;

				/*allow loss of biomass*/
				LB('sebiomasssi') = -M;
				UB('sebiomasssi') = 0;
				
			/*if has seed mass and seeds allow metabolic reactions*/
			ELSE
			
				/*allow flux in the forward direction for forward reactions*/
				UB(jseed)$(rxntype(jseed) = 1) = M;
				LB(jseed)$(rxntype(jseed) = 1) = 0;

				/*allow flux only in reverse direction for backwards reactions*/
				UB(jseed)$(rxntype(jseed) = -1) = 0;
				LB(jseed)$(rxntype(jseed) = -1) = -M;

				/*allow flux in both directions for reversible reactions*/
				UB(jseed)$(rxntype(jseed) = 0) = M;
				LB(jseed)$(rxntype(jseed) = 0) = -M;
			
			);
			
			omega_tm2 = omega_tm1;
			omega_tm1 = omega_t;
			omega_t = (sum(n1, curr_step(n1) * x_n(n1,'root')) * Y_thetat('leaf')) / (sum(n1, curr_step(n1) * x_n(n1,'leaf')) * Y_thetat('root'));
			
			eta_tm2 = eta_tm1;
			eta_tm1 = eta_t;
			eta_t = ((sum(n1, curr_step(n1) * x_n(n1,'seed')) * Y_thetat('leaf')) / (sum(n1, curr_step(n1) * x_n(n1,'leaf')) * Y_thetat('seed')))$(Y_thetat('seed') ne 0) + (0)$(Y_thetat('seed') eq 0);
			
			lambda_tm2 = lambda_tm1;
			lambda_tm1 = lambda_t;
			lambda_t = (sum(n1, curr_step(n1) * x_n(n1,'stem')) * Y_thetat('leaf')) / (sum(n1, curr_step(n1) * x_n(n1,'leaf')) * Y_thetat('stem'));
			
			lstarchrate = (A_l * sin(f_l * (t_n + b_l)))$(t_n > stage_0_7) + 0;
			ststarchrate = (A_st1 * sin(f_st1 * (t_n + b_st1)))$(t_n > stage_0_7) + 0;
			stsucroserate = (A_st2 * sin(f_st2 * (t_n + b_st2)))$(t_n > stage_0_7) + 0;
			
			/*LB('lbiomasssi') = 0;*/
			/*LB('rbiomasssi') = 0;*/
			/*LB('sebiomasssi') = 0;*/
			/*LB('stbiomasssi') = 0;*/
			
			/*solve the model to get the rates of growth*/
			SOLVE PATH780 USING LP MAXIMIZING obj;
			
			/*save the flux rates for this solution*/
			v_n(N,J) = v.l(J);
			
			/*calculate derivative estimate*/
			ds_dt = (3 * s_t - 4 * s_tm1 + s_tm2) / (2 * (c('2') - c('1')));
			dlnomegat_dt = (3 * log(omega_t) - 4 * log(omega_tm1) + log(omega_tm2)) / (2 * (c('2') - c('1')));
			dlnetat_dt = (3 * ((log(eta_t))$(eta_t ne 0) + 0) - 4 * ((log(eta_tm1))$(eta_tm1 ne 0) + 0) + ((log(eta_tm2))$(eta_tm2 ne 0) + 0)) / (2 * (c('2') - c('1')));
			dlnlambdat_dt = (3 * log(lambda_t) - 4 * log(lambda_tm1) + log(lambda_tm2)) / (2 * (c('2') - c('1')));
			
			/*calculate parameters for derivative estimates*/
			zeta_t = (c_theta('root') * (c_theta('leaf') * s_t + sum(n1, curr_step(n1) * x_n(n1,'leaf'))) - c_theta('leaf') * (c_theta('root') * s_t + sum(n1, curr_step(n1) * x_n(n1,'root')))) / (power(c_theta('leaf'),2) * power(s_t,2) + 2 * c_theta('leaf') * s_t * sum(n1, curr_step(n1) * x_n(n1,'leaf')) + sum(n1, curr_step(n1) * x_n(n1,'leaf')) * sum(n1, curr_step(n1) * x_n(n1,'leaf')));
			rho_t = (c_theta('seed') * (c_theta('leaf') * s_t + sum(n1, curr_step(n1) * x_n(n1,'leaf'))) - c_theta('leaf') * (c_theta('seed') * s_t + sum(n1, curr_step(n1) * x_n(n1,'seed')))) / (power(c_theta('leaf'),2) * power(s_t,2) + 2 * c_theta('leaf') * s_t * sum(n1, curr_step(n1) * x_n(n1,'leaf')) + sum(n1, curr_step(n1) * x_n(n1,'leaf')) * sum(n1, curr_step(n1) * x_n(n1,'leaf')));
			delta_t = (c_theta('stem') * (c_theta('leaf') * s_t + sum(n1, curr_step(n1) * x_n(n1,'leaf'))) - c_theta('leaf') * (c_theta('stem') * s_t + sum(n1, curr_step(n1) * x_n(n1,'stem')))) / (power(c_theta('leaf'),2) * power(s_t,2) + 2 * c_theta('leaf') * s_t * sum(n1, curr_step(n1) * x_n(n1,'leaf')) + sum(n1, curr_step(n1) * x_n(n1,'leaf')) * sum(n1, curr_step(n1) * x_n(n1,'leaf')));
			
			xi_t = sum(n1, curr_step(n1) * x_n(n1,'root')) * v.l('rbiomasssi') + sum(n1, curr_step(n1) * x_n(n1,'seed')) * v.l('sebiomasssi') + sum(n1, curr_step(n1) * x_n(n1,'stem')) * v.l('stbiomasssi') + sum(n1, curr_step(n1) * x_n(n1,'leaf')) * (v.l('rbiomasssi') * zeta_t + v.l('sebiomasssi') * rho_t + v.l('stbiomasssi') * delta_t) * ds_dt + sum(n1, curr_step(n1) * x_n(n1,'root')) * dlnomegat_dt +  sum(n1, curr_step(n1) * x_n(n1,'seed')) * dlnetat_dt +  sum(n1, curr_step(n1) * x_n(n1,'stem')) * dlnlambdat_dt; 
			k(N) = ((exp(v.l('lbiomasssi')) * sum(N1,Y_thetan(N1,'leaf'))) / sum(n1, curr_step(n1) * x_n(n1,'leaf'))) * (sum(n1, curr_step(n1) * x_n(n1,'leaf')) * v.l('lbiomasssi') + dmu_l_dt + xi_t);
			
			/*add the current element to the previous */
			Nm1(N) = YES;
			
			/*report results to troubleshooting file*/
			TROUBLESHOOT.ap = 1;
			PUT TROUBLESHOOT;
			TROUBLESHOOT.nr = 1;
			PUT /;
			PUT "********************************************************************************"/;
			PUT "*DAY: "Day.tl,system.tab,"TIME: ",hour.tl,system.tab,"STEP: ",N.tl,system.tab,system.tab,system.tab,system.tab,system.tab,"*"/;
			PUT "********************************************************************************"/;
			PUT "INPUT PARAMETERS"/;
			PUT "----------------"/;
			PUT "          t_n: ",t_n/;
			PUT "     past_0_7: ",past_0_7/;
			PUT "     past_1_2: ",past_1_2/;
			PUT "    diurnal_t: ",diurnal_t/;
			TROUBLESHOOT.nr = 2;
			PUT "       Y_n(n): ",Y_n(n)/;
			TROUBLESHOOT.nr = 1;
			PUT "         b(n): ",b(n)/;
			PUT "         c(n): ",c(n)/;
			PUT "          s_t: ",s_t:0:4,system.tab,"     s_tm1: ",s_tm1:0:4,system.tab,"     s_tm2: ",s_tm2:0:4/;
			TROUBLESHOOT.nr = 2;
			PUT "      omega_0: ",omega_t,system.tab," omega_tm1: ",omega_tm1,system.tab," omega_tm2: ",omega_tm2/;
			PUT "        eta_0: ",eta_t,system.tab,"   eta_tm1: ",eta_tm1,system.tab,"   eta_tm2: ",eta_tm2/;
			PUT "     lambda_0: ",lambda_t,system.tab,"lambda_tm1: ",lambda_tm1,system.tab,"lambda_tm2: ",lambda_tm2/;
			PUT "       x_leaf: ",x_n(N,'leaf')/;
			PUT "       x_root: ",x_n(N,'root')/;
			PUT "       x_seed: ",x_n(N,'seed')/;
			PUT "       x_stem: ",x_n(N,'stem')/;
			PUT "  lstarchrate: ",lstarchrate/;
			PUT " ststarchrate: ",ststarchrate/;
			PUT "stsucroserate: ",stsucroserate/;
			
			PUT /;
			PUT "SOLUTION"/;
			PUT "----------------------------------------------------------"/;
			PUT "                            OBJECTIVE (h^-1): ",obj.l/;
			PUT "                          LEAF GROWTH (h^-1): ",v.l('lbiomasssi')/;
			PUT "                          ROOT GROWTH (h^-1): ",v.l('rbiomasssi')/;
			PUT "                          SEED GROWTH (h^-1): ",v.l('sebiomasssi')/;
			PUT "                          STEM GROWTH (h^-1): ",v.l('stbiomasssi')/;
			TROUBLESHOOT.nr = 1;
			PUT "                                   MODELSTAT: ",PATH780.ModelStat/;
			PUT "                                   SOLVESTAT: ",PATH780.SolveStat/;
			TROUBLESHOOT.nr = 2;
			PUT "        PHOTOSYNTHESIS (LIGHT UPTAKE mmol/h): ", (v.l('R163ex[l]') * Y_thetat('leaf'))/;
			PUT "PHOTOSYNTHESIS (LIGHT UPTAKE mmol/gDWleaf*h): ", v.l('R163ex[l]')/;
			PUT "RESPIRATION (OX-PHOS RATE OF O2 USE, mmol/h): ", 0/;
			PUT "                TOTAL OXYGEN EXPORT (mmol/h): ", ((v.l('R044ex[r]') * Y_thetat('root')) + (v.l('R116ex[se]') * Y_thetat('seed')) + (v.l('R059ex[st]') * Y_thetat('stem')) + (v.l('R161ex[l]') * Y_thetat('leaf')))/;
			PUT "         LEAF OXYGEN EXPORT (mmol/gDWleaf*h): ", v.l('R161ex[l]'),system.tab,"LB: ",LB('R161ex[l]'),system.tab,"UB: ",UB('R161ex[l]')/;
			PUT "         ROOT OXYGEN EXPORT (mmol/gDWroot*h): ", v.l('R044ex[r]'),system.tab,"LB: ",LB('R044ex[r]'),system.tab,"UB: ",UB('R044ex[r]')/;
			PUT "         SEED OXYGEN EXPORT (mmol/gDWseed*h): ", v.l('R116ex[se]'),system.tab,"LB: ",LB('R116ex[se]'),system.tab,"UB: ",UB('R116ex[se]')/;
			PUT "         STEM OXYGEN EXPORT (mmol/gDWstem*h): ", v.l('R059ex[st]'),system.tab,"LB: ",LB('R059ex[st]'),system.tab,"UB: ",UB('R059ex[st]')/;
			PUT "                   TOTAL CO2 EXPORT (mmol/h): ", ((v.l('R164ex[l]') * Y_thetat('leaf')) + (v.l('R043ex[r]') * Y_thetat('root')) + (v.l('R125ex[se]') * Y_thetat('seed')) + (v.l('R066ex[st]') * Y_thetat('stem')))/;
			PUT "            LEAF CO2 EXPORT (mmol/gDWleaf*h): ", v.l('R164ex[l]'),system.tab,"LB: ",LB('R164ex[l]'),system.tab,"UB: ",UB('R164ex[l]')/;
			PUT "            ROOT CO2 EXPORT (mmol/gDWroot*h): ", v.l('R043ex[r]'),system.tab,"LB: ",LB('R043ex[r]'),system.tab,"UB: ",UB('R043ex[r]')/;
			PUT "            SEED CO2 EXPORT (mmol/gDWseed*h): ", v.l('R125ex[se]'),system.tab,"LB: ",LB('R125ex[se]'),system.tab,"UB: ",UB('R125ex[se]')/;
			PUT "            STEM CO2 EXPORT (mmol/gDWstem*h): ", v.l('R066ex[st]'),system.tab,"LB: ",LB('R066ex[st]'),system.tab,"UB: ",UB('R066ex[st]')/;
			PUT "                                  k_m (h^-1): ", k_m:0:8/;
			PUT "                           v_k_m leaf (h^-1): ", v.l('lbiomass_main'):0:8/;
			PUT "                                  k_s (h^-1): ", k_s:0:8/;
			PUT "                           v_k_m leaf (h^-1): ", v.l('lbiomass_scen'):0:8/;
			
			PUT /;
			PUT "POST SOLUTION PARAMETERS"/;
			PUT "---------------------------"/;
			PUT "        ds_dt: ",ds_dt/;
			PUT " dlnomegat_dt: ",dlnomegat_dt/;
			PUT "   dlnetat_dt: ",dlnetat_dt/;
			PUT "dlnlambdat_dt: ",dlnlambdat_dt/;
			PUT "       zeta_t: ",zeta_t/;
			PUT "        rho_t: ",rho_t/;
			PUT "      delta_t: ",delta_t/;
			PUT "         xi_t: ",xi_t/;
			PUT "          k_n: ",k(n)/;
			
			IF((PATH780.ModelStat eq 4),
			
				/*kill the run now*/
				/*abort if there is an infeasible point because anything after that point will be nonsense*/
				ABORT "ABORTED DUE TO INFEASIBLE MODEL STATUS";
			
			);
			
			/*for certain time points do FVA quick to make sure everything is working right*/
			IF(((t_n eq 1) OR (t_n eq 70) OR (t_n eq 90) OR (t_n eq 90) OR (t_n eq 177) OR (t_n eq 181) OR (t_n eq 770) OR (t_n eq 810) OR (t_n eq 1155) OR (t_n eq 1170) OR (t_n eq 1190) OR (t_n eq 1199)),
			
				PUT //;
				
				PUT "FVA RESULTS FOR GIVEN TIMEPOINT"/;
				PUT "REACTION",system.tab,system.tab,system.tab,system.tab,system.tab,"LB",system.tab,"   MINIMUM",system.tab,"   MAXIMUM",system.tab,system.tab,system.tab,"UB"/;
				
				PUTCLOSE;
				
				LOOP(J1,
				
					/*fix biomass rates to see variablility at the current growth rate*/
					v.fx('lbiomasssi') = v.l('lbiomasssi');
					v.fx('rbiomasssi') = v.l('rbiomasssi');
					v.fx('sebiomasssi') = v.l('sebiomasssi');
					v.fx('stbiomasssi') = v.l('stbiomasssi');
				
					TROUBLESHOOT.ap = 1;
					TROUBLESHOOT.lw = 18;
					PUT TROUBLESHOOT;
				
					sigma(J) = 0;
					sigma(J1) = 1;
					PUT J1.tl;
					
					TROUBLESHOOT.lw = 12;
					
					PUT LB(J1);
					
					/*SOLVE FVA_PATH780 USING LP MINIMIZING obj;
					PUT obj.l;
					
					SOLVE FVA_PATH780 USING LP MAXIMIZING obj;
					PUT obj.l;*/
				
					PUT UB(J1);
					
					PUT /;
					
					PUTCLOSE;
					
					/*unfix biomass*/
					v.lo('lbiomasssi') = -M;
					v.up('lbiomasssi') = M;
					v.lo('rbiomasssi') = -M;
					v.up('rbiomasssi') = M;
					v.lo('sebiomasssi') = -M;
					v.up('sebiomasssi') = M;
					v.lo('stbiomasssi') = -M;
					v.up('stbiomasssi') = M;
				
				);
			
				/*ABORT "ABORT FOR TROUBLESHOOTING SOLUTION POINT";*/
			
			);
			
		);
		
		/*calculate the mass step*/
		dY_dt = 1 * sum(N1, b(N1) * k(N1));
		
		/*write mass sheet line*/
		MASSSHEET.ap = 1;
		PUT MASSSHEET;
		MASSSHEET.pc = 5;
		MASSSHEET.nr = 1;
		
		PUT day.tl,hour.tl,t_0;
		
		MASSSHEET.nr = 2;
		PUT Y_0,Y_thetat('leaf'),Y_thetat('root'),Y_thetat('seed'),Y_thetat('stem');
		
		LOOP(N,
		
			PUT Y_n(N);
		
		);
		
		PUT dY_dt/;
		
		MASSSHEET.nr = 1;
		PUTCLOSE;
		
		Y_Deltat = Y_0 + Y_G * dY_dt;
		Gamma_t(j) = (c('2') - c('1')) * 1 * (2 * sum(N, v_n(N,J)) - v_n('1',J) + ((sum(N, v_n(N,J))) / (sum(N, 1)))) / (2 * sum(N, 1)); 
		Z_Deltat(I) = Z_0(I) + sum(J, S(I,J) * Gamma_t(J));
		
		TROUBLESHOOT.ap = 1;
		PUT TROUBLESHOOT;
		PUT /;
		PUT "SOLUTION SUMMARY"/;
		PUT "----------------------------"/;
		TROUBLESHOOT.nr = 1;
		PUT "CURRENT TIME: ",t_0/;
		PUT "   NEXT TIME: ",t_Deltat/;
		TROUBLESHOOT.nr = 2;
		PUT "dY/dt_rk est: ", dY_dt/;
		PUT "   NEXT MASS: ",Y_Deltat/;
		PUT /;
		
		PUTCLOSE;
		
		/*update for t_0 and Y_0 for the next step*/
		t_0 = t_Deltat;
		Y_0 = Y_Deltat;
		Z_0(I) = Z_Deltat(I);
		
		/*write flux output file*/
		FLUXSHEET.ap = 1;
		PUT FLUXSHEET;
		PUT Day.tl,hour.tl,t_0;
		
		LOOP(J,
		
			PUT v_n('1',J);
		
		);
		
		PUT /;
		
		PUTCLOSE;
		
		/*If a negative mass or concentration is returned this is due to the slack on the*/
		/*mass balance constraint these values are generally less than epsilon so will*/
		/*reset them to zero whe reporting them*/
		LOOP(I,
		
			IF((Z_0(I) < 0),
			
				Z_0(I) = 0;
			
			);
		
		);
		
		CONCSHEET.ap = 1;
		PUT CONCSHEET;
		PUT Day.tl,hour.tl,t_0;
		
		LOOP(I,
		
			PUT Z_0(I);
		
		);
		
		PUT /;
		
		IF ((Y_Deltat < 0),
		
			TROUBLESHOOT.ap = 1;
			PUT TROUBLESHOOT;
			ABORT "Negative plant mass returned";
			PUTCLOSE;
			
		);
		
		/*if max time has been reached kill the run*/
		IF ((t_n ge Max_Time),
		
			stop = 1;
		
		);
		
	);
	
	/*if max time has been reached kill the run*/
	IF ((t_n ge Max_Time),
		
		stop = 1;
		
	);
	
);

TROUBLESHOOT.ap = 1;
PUT TROUBLESHOOT;

PUT "SUCCESSFUL COMPLETION!";

PUTCLOSE;
