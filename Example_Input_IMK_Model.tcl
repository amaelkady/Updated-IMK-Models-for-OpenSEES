# ####################################################################################
# Sample Model of a ZeroLength Spring with an IMKPeakOriented Response
#
# Ahmed Elkady, 05-05-2020
# #####################################################################################
wipe all

# Define model
model BasicBuilder -ndm 2 -ndf 3

# 2 Nodes are needed to define a rotational spring
node 1 0.0 0.0 0.0
node 2 0.0 0.0 0.0

# Fixed Conditions for one of the two nodes
fix 1 1 1 1
fix 2 1 0 1 
# #####################################################################################
# #####################################################################################
# #####################################################################################
# #####################################################################################

# IMPUT
set Ke 		2.e7;
set My 		2.e4;
set Mc_My 	1.5;
set thetap  0.03;
set thetapc 0.10;
set thetau  0.15;
set lambda  1.0;
set Mres_My 0.1;


uniaxialMaterial IMKPeakOriented 1 $Ke $thetap $thetapc $thetau $My $Mc_My $Mres_My $thetap $thetapc $thetau $My $Mc_My $Mres_My $lambda $lambda $lambda $lambda 1 1 1 1 1 1;

# #####################################################################################
# #####################################################################################
# #####################################################################################
# #####################################################################################
# #####################################################################################

# Disp History
set disp [list 0.0 0	0.375	-0.375	0.5	-0.5 0.75	-0.75	1	-1	1.5	-1.5	2	-2	3	-3	4	-4];

# Define Zero Length Rotational Spring
element zeroLength  1 1 2 -mat 1 -dir 2
	
# Unit Load Pattern in Rotational Degree of Freedom
pattern Plain 1 "Linear" {
		# nd    FX  FY
	load 2     0.0 1.0 0.0
}

# Recorders # ------------------------------------------------------------------------
recorder Element -file Force.out -ele 1 force
recorder Element -file Disp.out  -ele 1 deformation

# Set analysis parameters # ----------------------------------------------------------
test EnergyIncr 1.0e-8 50 0
algorithm Newton
system UmfPack
numberer RCM
constraints Plain

set scale 0.03;  # scale factor for the disp history
set NSteps 100; # controls resolution of disp increment 

set LoopLength [llength $disp]
set h 1
set controlNode 2

# Run the static cyclic analysis
set dU1 0;
set controlNodeDisp [nodeDisp $controlNode 2]

while {$h < $LoopLength} {
	# List is zero index
	set index [expr $h]
	# We need to relative deformation of the loading protocol
	# Substract dUi+1 - dUi
	set D1 [lindex $disp $index];
	set D2 [lindex $disp $index-1];	
	set dU1 [expr -($D1-$D2)*$scale];
	
	# Create Nsteps from Amplitude to Amplitude
	set dU [expr ($dU1)/$NSteps]
	
	# Displacement Control Integrator
	integrator DisplacementControl $controlNode 2 $dU 1 $dU $dU
	analysis Static

	set ok [ analyze $NSteps]
	set h [expr $h + 1 ]
}

wipe all;