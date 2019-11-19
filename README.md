# PlaneMaker
Read Me File for Plane Maker 2

Plane Maker calculates performance aspects of an aircraft using basic design parameters. The user can input information in lines 29-85.
From there the code will process the data and perform the necessary calculations. Then, outputs are delivered in the terminal. This program
has a special design subroutine. This subroutine takes the data you've entered and uses it to come up with a very accurate weight estimate.
It does this using a empty weight fraction that has been determined from historical trends. The subroutine outputs recomended weights,
airspeeds, and ensures landing location is ok. This subroutine is toggleable using a boolean in a design mode controls input section.

Steps:
1) Find file for airfoil information
2) Process file for program reading (explained below)
3) Open program
4) Enter all the information for your aircraft in imperial units
5) To use the Design Mode, change the boolean next to 'RunSysDesign' to True
6)Run the program
7) Enter the values given by the Design Mode
8) Rerun Design Mode to ensure accuracy
9) Change the boolean next to 'RunSysDesign' to false
10) Run the program
11) Evaluate the outputs and make changes to your aircraft as needed.
12) If you encounter any issues, report them to gjdrangstveit@gmail.com with subject 
'Plane Maker 2 Issue'

Airfoil File:
When importing a file for containing airfoil data, it must follow a certain format. All values should be arranged in columns, with angle
of attack on the left and lift coefficients on the right seperated by commas. Example:
-9.75,-0.5187
-9.5,-0.5473
-9.25,-0.5826
-9,-0.6085
-8.75,-0.604
-8.5,-0.5974
When inputting the file location into the program, be sure to change an instances of '\' to '/'. If this step is not done the code will not
be able to access the file. Example:
C:\Users\mainuser\Documents\PythonStuff\Airfoil.txt
becomes
C:/Users/mainuser/Documents/PythonStuff/Airfoil.txt

Empty Weight Ratio:
Below is a list of ratios taken from Aircraft Design: A Conceptual Approach.

Sailplane-unpowered:              0.590-0.640
Sailplane-powered:                0.575-0.660
Homebuilt-metal/wood:             0.625-0.660
Homebuilt-composite:              0.525-0.575
General Aviation-single engine:   0.525-0.620
General Aviation-twin engine:     0.600-0.651
Agricultural aircraft:            0.551-0.580
Twin turboprob:                   0.550-0.610
Flying boat:                      0.690-0.730
Jet trainer:                      0.620-0.660
Jet fighter:                      0.500-0.650
Military cargo/bomber:            0.350-0.445
Jet transport:                    0.450-0.555

Use a ratio within the range given for the type of aircraft you are designing. If you expect your aircraft to be lighter, use the smaller
value or the larger value for larger weights.
Input this value for 'EmptyWeightRatio'.

Lift Coefficient With Flaps:
Below is the list of possible lift coefficients for landing and take off. Multiply this value by the cosine of your sweep angle and enter
that as your input.

High Lift Device                   Takeoff   Landing
Plain flap                         1.4-1.6   1.7-2.0
Single-slotted flap                1.5-1.7   1.8-2.2
single-slotted fowler              2.0-2.2   2.5-2.9
double-slotted fowler              1.7-1.95  2.3-2.7
double-slotted fowler & slats      2.3-2.6   2.8-3.2
triple-slotted fowler & slats      2.4-2.7   3.2-3.5

If you aren't sure what value to use, use the lowest value in the range.

Weight Estimation Explanation:
The program generates an accurate estimate of your weight. To do this, it uses historical data about the weights of planes that was
used to generate ranges for the empty weight fraction of aircraft. The program first runs through a loop and calculates all fuel used during
the flight. Then it finds the total weight of the craft, including any user input weights.  It then multiplies the weight of the aircraft
by a historical ratio to get an estimate of the weight of the craft. Then, it runs through a loop again and once again calculates the sum
of the weights. It does this repeatedly until the estimate that it generates in the loop is within a given accuracy to the value the loop
used.
