# START OF PROGRAM



#imports math library and plotting library
import math
import matplotlib.pyplot as plt



# PROGRAM INFORMATION



#PLANE MAKER 2
#CREATED BY GUNNAR DRANGSTVEIT
#v<2.4.2>

#Read the ReadMe file that came with this program for instructions and related information.
#Design Mode: Change 'RunSysDesign' to 'True' if you wish to use Design Mode, and 'False' if you do not.
#Output Modes: Change 'AbridgedOutputs' to 'True' if you want shorter outputs, and 'False' if you do not.


# USER INPUTS (ENTER ALL IN DESIGNATED IMPERIAL UNITS):



#fuselage dimensions
FuselageLength = 123 #ft
FuselageDiameter = 123 #ft

#wing dimensions
RootChord = 123 #ft
TipChord = 123 #ft
WingSpan = 123 #ft
SweepAngle = 123 #deg
LocationOfMaxThickness = 0.123 #unitless
MaxThickness = 0.123 #unitless

#file containing lift coifficients and AoA's
FileName = 'C:/Users/drago/Documents/Python Scripts/BACXXX_RE1milN9.txt' #Change instances of '\' to '/'

#cruise phase informations
CruiseAltitude = 12345 #ft
CruiseAirspeed = 123 #fps
FlightRange = 1234 #mi

#aircraft weight information
WeightOfFuselage = 123456 #lb
WeightOfPerson = 170 #lb
WeightOfCargo = 12345 #lb
WeightOfFuel = 123456 #lb
WeightOfFuelReserve = 12345 #lb
WeightOfEngines = 1234 #lb

#onboard quantities
People = 123 #people
Engines = 1 #unitless

#jet engine information
EngineThrust = 12345 #lb (per engine)
SpecificFuelConsumption = 0.123 #1/hr
ReverseThrust = 12345 #lb
EngineLength = 12 #ft
EngineDiameter = 12 #ft

#information of takeoff, climb, and landing locations
TakeoffAltitude = 0 #ft
TakeoffRunwayLength = 12345 #ft
ClimbStep = 1000 #ft (sets accuracy of climb sim)
UnpoweredTestAltitude = 12345 #ft
LandingAltitude = 0 #ft
LandingRunwayLength = 12345 #ft
PilotReactionTime = 3 #sec
LandingLiftCoefficentWithFlaps = 1.23 #unitless
TakeoffLiftCoefficientWithFlaps = 1.23 #unitless

#iterative Design Mode controls
RunSysDesign = False #boolean (True or False)
WeightAccuracy = 0.001 #lb (lower than 0.001 may take long run times)
EmptyWeightRatio = 0.123 #unitless (Found in Conceptual Design Textbook)

#output controls
AbridgedOutputs = False #boolean (True or False)



# END OF USER INPUTS



#/////////////////////////////////////////////////// WARNING: CHANGING ANY OF THE FOLLOWING CODE CAN AND WILL RESULT IN INCORRECT RESULTS, DO NOT TAMPER \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#/////////////////////////////////////////////////// WARNING: CHANGING ANY OF THE FOLLOWING CODE CAN AND WILL RESULT IN INCORRECT RESULTS, DO NOT TAMPER \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#/////////////////////////////////////////////////// WARNING: CHANGING ANY OF THE FOLLOWING CODE CAN AND WILL RESULT IN INCORRECT RESULTS, DO NOT TAMPER \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#/////////////////////////////////////////////////// WARNING: CHANGING ANY OF THE FOLLOWING CODE CAN AND WILL RESULT IN INCORRECT RESULTS, DO NOT TAMPER \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#/////////////////////////////////////////////////// WARNING: CHANGING ANY OF THE FOLLOWING CODE CAN AND WILL RESULT IN INCORRECT RESULTS, DO NOT TAMPER \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#/////////////////////////////////////////////////// WARNING: CHANGING ANY OF THE FOLLOWING CODE CAN AND WILL RESULT IN INCORRECT RESULTS, DO NOT TAMPER \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#/////////////////////////////////////////////////// WARNING: CHANGING ANY OF THE FOLLOWING CODE CAN AND WILL RESULT IN INCORRECT RESULTS, DO NOT TAMPER \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#



# START OF FILE PARSING AND AERODYNAMIC CALCULATIONS



#opens the file and reads it into list
file = open(FileName, 'r')
airFoil = file.read()
#splits the list of values at each line break
newFoil = airFoil.split('\n')

#sets up values needed for loop
i=0
clpts = []
apts = []
#creates list of lists for each data point, converts strings to floats
for line in newFoil:
    #breaks each line at the seperating space
    nline = line.split(',')
    #creates new sublist in list of data points
    clpts = clpts + [float(nline[1])]
    apts = apts + [float(nline[0])]
    #increments counter
    i+=1

#sets counter value for loop
j=0
#finds where alpha is zero
for num in apts:
    if num == 0:
        #records corresponding coefficient
        cl0 = clpts[j]
        break
    #increments counter
    j+=1

#sets values needed for loop
k = 0
small = 100
#finds where lift coefficient is closest to zero
for num in clpts:
    if abs(num) < abs(small):
        #sets new global minimum
        small = num
        #records corresponding alpha
        a0 = apts[k]
    #increments counter
    k+=1
#sets lift coefficient to minimum value found by the loop
cl = small

#total weight of aircraft
Weight = WeightOfFuselage + WeightOfFuel + WeightOfFuelReserve + (WeightOfPerson*People) + WeightOfCargo + (WeightOfEngines*Engines)
#average chord of wing
chord = (RootChord+TipChord)/2
#slope of airfoil lift curve
a_0 = (cl0-cl)/(0-a0)
#aspect ratio
AR = WingSpan/chord

#calculates oswald efficiency factor
if SweepAngle == 0:
    #for unswept wing
    OEF = 1.78*(1-0.045*pow(AR,0.68))-0.64
elif SweepAngle > 30:
    #for sweep angle over 30deg
    OEF = 4.61*(1-0.045*pow(AR,0.68))*pow(math.cos(math.radians(SweepAngle)),0.15)-3.1
elif SweepAngle <= 30 and SweepAngle > 0:
    #for sweep angle between 0-30deg
    OEF = ((30-SweepAngle)/30)*(1.78*(1-0.045*pow(AR,0.68))-0.64) + (SweepAngle/30)*((4.61*(1-0.045*pow(AR,0.68))*pow(math.cos(math.radians(SweepAngle)),0.15))-3.1)
#induced drag factor
K = 1/(3.14159*OEF*AR)
#slope of the wing lift curve
a = (a_0*math.cos(math.radians(SweepAngle)))/(1+K*a_0*math.cos(math.radians(SweepAngle)))

#sets up empty list for loop
CLpts = []
#calculates values of wing lift coefficents
for x in range(10):
    CLpts = CLpts + [a*x + cl0]

#function that finds temperature at a given altitude
#found in kelvin then converted because I couldn't be bothered to change it (legacy code)
def GetTemperature(alt):
    #converts height in feet to kilometers
    height = alt/3280.84
    if height >=0 and height <=11:
        #for 0-11km
        temp = 288.16 + (-6.5*(height-0))
    elif height > 11 and height <= 25:
        #for 11-25km
        temp = 288.16 + (-6.5*(11-0))
    elif height > 25 and height <= 47:
        #for 25-47km
        temp = (288.16+(-6.5*(11-0)))+(3*(height-25))
    #converts temperature in Kelvin to Rankine
    temp = 1.8*temp
    return temp

#function that calculates density from a given altitude
#found in kg/m^3 then converted to slug/ft^3 becuase I didn't want to change it (legacy code)
def GetDensity(alti):
    #gets temperature and converts to Kelvin
    temp = GetTemperature(alti)*(5/9)
    #converts altitude in feet to kilometers
    nalt = alti/3280.84
    if nalt >= 0 and nalt <= 11:
        #for 0-11km
        dense = 1.225*(pow((temp/288.16),((-9.80065/((-0.0065)*287))-1)))
    elif nalt > 11 and nalt <= 25:
        #for 11-25km
        ntemp = 288.16+((-6.5)*(11-0))
        ndense = 1.225*(pow((ntemp/288.16),((-9.80065/((-0.0065)*287))-1)))
        dense = ndense*(math.exp((-9.80065/(287*ntemp))*((alti*1000) - 11000)))
    elif nalt > 25 and nalt <= 45:
        #for 25-47km
        ntemp = 288.16+((-6.5)*(11-0))
        ntemp2 = ntemp + ((-6.5)*(11-0))
        ndense = 1.225*(pow((ntemp/288.16),((-9.80065/((-0.0065)*287))-1)))
        ndense2 = ndense*(math.exp((-9.80065/(287*ntemp))*((25*1000)-11000)))
        dense = ndense2*(pow((ntemp2/ntemp),((-9.80065/((0.003)*287))-1)))
    #converts density from kg/m^3 to slug/ft^3
    Rdense = dense/515.379
    return Rdense

#Reynolds Number of wing and fuselage
Rewing = (GetDensity(CruiseAltitude)*CruiseAirspeed*chord)/0.000000374
Refuselage = (GetDensity(CruiseAltitude)*CruiseAirspeed*FuselageLength)/0.000000374
Renacelle = (GetDensity(CruiseAltitude)*CruiseAirspeed*EngineLength)/0.000000374
#speed of sound at cruise altitude
SoS = math.sqrt(1.4*1716*GetTemperature(CruiseAltitude))
#mach number
mach = CruiseAirspeed/SoS
#coefficient of friction for wing
lamflow = 1.328/math.sqrt(Rewing)
turbflow = (0.455/pow(math.log10(Rewing),2.58))*(1/pow((1+0.144*pow(mach,2)),0.65))
cfwing = (0.1*lamflow)+(0.9*turbflow)
#coefficient of friction for fuselage
cffuselage  = (0.455/pow(math.log10(Refuselage),2.58))*(1/pow((1+0.144*pow(mach,2)),0.65))
#coefficient of friction for the engine nacelle
cfnacelle  = (0.455/pow(math.log10(Renacelle),2.58))*(1/pow((1+0.144*pow(mach,2)),0.65))
#form factor for wing
FFwing = (1+(0.6/LocationOfMaxThickness)*(MaxThickness)+100*pow(MaxThickness,4))*(1.34*pow(mach,0.18)*pow(math.cos(math.radians(SweepAngle)),0.28))
#form factor of fuselage
ffuselage = FuselageLength/FuselageDiameter
FFfuselage = 0.9 + (5/pow(ffuselage,1.5))+(ffuselage/400)
#form factor for engine nacelle
fnacelle = EngineLength/math.sqrt((4/3.14159)*(3.14159*pow(FuselageDiameter/2,2)))
FFnacelle = 1+(0.35/fnacelle)
#wing area
S = WingSpan*chord
#surface area of fuselage
SAcone = 3.14159*(FuselageDiameter/2)*((FuselageDiameter/2)+math.sqrt(pow((FuselageLength*0.1),2)+pow(FuselageDiameter/2,2)))
SAcylinder = 2*3.14159*(FuselageDiameter/2)*(FuselageLength*0.9)
Swet = SAcone + SAcylinder
#surface area of nacelle
SAnacelle = (2*3.14159*(EngineDiameter/2)*(EngineLength))+(2*3.14159*pow(EngineDiameter/2,2))
#no lift drag
cd0wing = (cfwing*FFwing*1.5*2)
cd0fuselage = (cffuselage*FFfuselage*1.5*(Swet/S))
cd0nacelle = (cfnacelle*FFnacelle*1.5*(SAnacelle/S))
cd0 = cd0wing + cd0fuselage + cd0nacelle

#creates values for loop
m = 0
big = 0
#finds largest airfoil lift coefficient and corresponding angle
for num in clpts:
    #checks if number is larger than global max
    if abs(num) > abs(big):
        #records new global max
        big = num
        #records corresponding angle
        amax = apts[m]
    #increments counter
    m += 1

#largest wing lift coefficient
CLmax = a*amax + cl0

#maximum aerodynamic efficiency
Emax = 1/(2*math.sqrt(K*cd0))

ndiff = 100
#finds wing mounting angle
for num in apts:
    clmount = math.sqrt(cd0/K)
    diff = abs(clmount-(a*num+cl0))
    if diff < ndiff:
        ndiff = diff
        mountang = num

#lift coefficient for max aerodynamic efficiency
CLEM = math.sqrt(cd0/K)
#gets corresponding angle of attack
ndiff = 100
#finds wing mounting angle
for num in apts:
    diff = abs(CLEM-(a*num+cl0))
    if diff < ndiff:
        ndiff = diff
        aEM = num
        


# END OF FILE PARSING AND AERODYNAMIC CALCULATIONS


# START OF DESIGN MODE



if RunSysDesign == True:
    #sets up iteration counter
    count = 1
    #sets up initial values for loop
    currWdiff = 100
    Weightempty = WeightOfFuselage
    Weightfuel1 = WeightOfFuel
    Weightfuelreserve1 = WeightOfFuelReserve
    velocity_sysdesign = CruiseAirspeed
    #loop runs while sum of inaccuracy is larger than user specified
    while currWdiff > WeightAccuracy:
        #finds total weight of aircraft
        Weight1 = Weightempty + Weightfuel1 + Weightfuelreserve1 + (WeightOfPerson*People) + WeightOfCargo + (WeightOfEngines*Engines)
        #creates values for loop
        alt = TakeoffAltitude
        tclimb = 0
        Wclimb = 0
        newWeight = Weight1
        #loop that calculates flight from take off location to cruise altitude
        while alt <= CruiseAltitude:
            #finds thrust at altitudes
            Talt1 = (GetDensity(alt)/GetDensity(0))*(EngineThrust*Engines)
            Talt2 = (GetDensity(alt+ClimbStep)/GetDensity(0))*(EngineThrust*Engines)
            #values that makes later calculations easier
            gamma1 = 1 + math.sqrt(1+((12*K*cd0)/pow((Talt1/newWeight),2)))
            gamma2 = 1 + math.sqrt(1+((12*K*cd0)/pow((Talt2/newWeight),2)))
            #fastest climb velocity
            Vfastestclimb1 = math.sqrt((((Talt1)/S)*gamma1)/(3*GetDensity(TakeoffAltitude)*cd0))
            Vfastestclimb2 = math.sqrt((((Talt2)/S)*gamma2)/(3*GetDensity(TakeoffAltitude)*cd0))
            #sine of the flight path angle
            singammamax1 = ((Talt1)/newWeight)*(1-(gamma1/6))-((6*K*cd0)/(gamma1*((Talt1)/newWeight)))
            singammamax2 = ((Talt2)/newWeight)*(1-(gamma2/6))-((6*K*cd0)/(gamma2*((Talt2)/newWeight)))
            #vertical velocity
            RoCmax1 = Vfastestclimb1*singammamax1
            RoCmax2 = Vfastestclimb2*singammamax2
            RoCavg = (RoCmax1+RoCmax2)/2
            #average thrust
            Tavg = (Talt1+Talt2)/2
            #time to climb
            dt = ClimbStep/RoCavg
            tclimb+=dt
            #fuel burned during climb
            dW = (SpecificFuelConsumption*(1/3600))*(Tavg*dt)
            #weight at end of step
            newWeight-=dW
            #weight loss during climb
            Wclimb+=dW
            #increments altitude
            alt+=ClimbStep
        #aircraft weight at start of cruise
        Wstartcruise = Weight1 - Wclimb

        
        #time of flight assuming constant altitude and constant velocity
        t = (FlightRange*5280)/velocity_sysdesign
        #estimated weight loss
        dW = SpecificFuelConsumption*(1/3600)*(EngineThrust*Engines*(GetDensity(CruiseAltitude)/GetDensity(0)))*t
        Wendcruise = Wstartcruise-dW
        #averageweight
        Wavg = (Wstartcruise+Wendcruise)/2
        
        #velocity for max efficiency
        VEM = math.sqrt(((2*(Wavg/S))/GetDensity(CruiseAltitude))*math.sqrt(K/cd0))
        #stall velocity for cruise
        Vstallcruise = math.sqrt((2*(Wavg/S))/(GetDensity(CruiseAltitude)*CLmax))
        #coffin velocity for cruise
        Vstarcruise = math.sqrt((2*(Wavg/S)*8)/(GetDensity(CruiseAltitude)*CLmax))
        #speed of sound at cruise altitude
        SoScruise = math.sqrt(1.4*1716*GetTemperature(CruiseAltitude))  
        #80% of speed of sound
        V80 = 0.8*SoScruise
        #compares velocity for max efficiency and makes sure it is possible
        if VEM > Vstallcruise and VEM < Vstarcruise:
            vsys = 1
        elif V80 > Vstallcruise and V80 < Vstarcruise:
            vsys = 1
        else:
            vsys = 2        
        
        
        # Start Weight Design Iteration Portion
        
        
        #weight of all parts of aircraft
        We1 = Weightempty
        Wf = Weightfuel1+Weightfuelreserve1
        Wp = WeightOfCargo
        Wc = WeightOfPerson*People
        Wn = WeightOfEngines*Engines
        #total aircraft weight
        W0 = We1+Wf+Wp+Wc+Wn
        #estimate of empty weight using total weight
        We = EmptyWeightRatio*W0
        #compares estimate of empty weight to actual value
        Wediff = Weightempty-We
        #total fuel loss during flight
        Wlosstotal = Wclimb + dW
        #reserve fuel required
        Wreserve = SpecificFuelConsumption*(EngineThrust*Engines*(GetDensity(15000)/GetDensity(0)))*0.75
        #compares flight fuel loss to fuel available
        Weightdiff = Weightfuel1 - Wlosstotal
        #compares required fuel reserve to provided
        Weightreservediff = Weightfuelreserve1 - Wreserve
        
        Weightempty = We
        Weightfuel1 = Wlosstotal
        Weightfuelreserve1 = Wreserve
        currWdiff = abs(Wediff)+abs(Weightdiff)
        
        
        # Start Thrust Recommendation Portion
        
        
        #required thrust at altitude cruise
        Treqcruisestart = (0.5*GetDensity(CruiseAltitude)*pow(velocity_sysdesign,2)*cd0*S)+((K*pow(Wstartcruise/S,2)*S)/(0.5*GetDensity(CruiseAltitude)*pow(velocity_sysdesign,2)))
        Treqcruiseend = (0.5*GetDensity(CruiseAltitude)*pow(velocity_sysdesign,2)*cd0*S)+((K*pow(Wendcruise/S,2)*S)/(0.5*GetDensity(CruiseAltitude)*pow(velocity_sysdesign,2)))
        #average weight during cruise
        Wavg = (Wstartcruise + Wendcruise)/2
        #velocity of minimum turn radius
        Vminradius = math.sqrt((4*K*(Wavg/S))/(GetDensity(CruiseAltitude)*((EngineThrust*Engines*(GetDensity(CruiseAltitude)/GetDensity(0)))/Wavg)))
        #load factor for minimum turn radius
        nminradius2 = 2-((4*K*cd0)/pow((EngineThrust*Engines*(GetDensity(CruiseAltitude)/GetDensity(0)))/Wavg,2))     
        #minimum Thrust for turning
        Tminturn = (0.5*GetDensity(CruiseAltitude)*pow(Vminradius,2)*cd0*S)+((K*pow(Wavg/S,2)*S*nminradius2)/(0.5*GetDensity(CruiseAltitude)*pow(Vminradius,2)))
        #finds largest thrust value
        if Tminturn >= Treqcruisestart and Tminturn >= Treqcruiseend:
            Tminimum = Tminturn
        elif Treqcruisestart >= Tminturn and Treqcruisestart >= Treqcruiseend:
            Tminimum = Treqcruisestart
        elif Treqcruiseend >= Tminturn and Treqcruiseend >= Treqcruisestart:
            Tminimum = Treqcruiseend
        
        
        # Start Runway Length Checker
        
        
        #finds max wing lift coefficient at takeoff
        CLmaxland = TakeoffLiftCoefficientWithFlaps*math.cos(math.radians(SweepAngle))
        #finds stall velocity at takeoff location
        Vstall = math.sqrt((2*(Wendcruise/S))/(GetDensity(TakeoffAltitude)*CLmaxland))        
        #calculates take off velocity
        Vtakeoff = 1.2*Vstall
        #takeoff acceleration
        atakeoff = 32.17*((EngineThrust*Engines)/Weight)
        #take off distance
        dtakeoff = pow(Vtakeoff,2)/(2*atakeoff)        
        
        #finds max wing lift coefficient at landing
        CLmaxland = LandingLiftCoefficentWithFlaps*math.cos(math.radians(SweepAngle))
        #finds stall velocity at landing location
        Vstall = math.sqrt((2*(Wendcruise/S))/(GetDensity(LandingAltitude)*CLmaxland))
        #touchdown velocity
        Vtouchdown = 1.15*Vstall
        #distance traveled while rolling
        sfreeroll = Vtouchdown*PilotReactionTime
        #some things that make later calculations easier
        JT = (ReverseThrust*Engines)/Wendcruise
        JA = (GetDensity(LandingAltitude)*cd0)/(2*(Wendcruise/S))
        #finds distance required to stop
        s = sfreeroll + (1/(2*32.2*JA))*math.log(1+((JA/JT)*pow(Vtouchdown,2)),2.71828)
        
        #increments counter
        count+=1
        
        #if statement for excessive run time error
        if count == 100:
            print('\nsysalarm: Run Count in Excess of 100 Iterations')
        elif count == 1000:
            print('\nsysalarm: Run Count in Excess of 1,000 Iterations')
        elif count == 10000:
            print('\nsysalarm: Run Count in Excess of 10,000 Iterations')
            print('sysalarm: Run Count Excessive, Self-Termination in Progress')
            exit()
    
    
    #outputs for Design System
    
    
    #number of loops done
    print('\nsysinfo:',count,'Iterations')
    #calculated fuselage weight
    print('\nsysdesign: Use Fuselage Weight',round(Weightempty,2),'lb')
    #calculated fuel weight
    print('\nsysdesign: Use Fuel Weight',round(Weightfuel1,2),'lb')
    #calculated fuel reserve weight
    print('\nsysdesign: Use Fuel Reserve',round(Weightfuelreserve1,2),'lb')
    #calculated minimum engine thrust
    print('\nsysdesign: Minimum Thrust Per Engine',round(Tminimum/Engines,2),'lb')
    #best velocity for cruise
    if vsys == 1:
        print('\nsysdesign: Use Cruise Speed',round(velocity_sysdesign,2),'fps')
    elif vsys == 2:
        print('\nsysdesign: Unable To Recommend Cruise Speed')
    else:
        print('\nsyserror: Error in Cruise Speed Calculations')
    #checks that runways are long enough for calculated flight
    if dtakeoff > TakeoffRunwayLength:
        print('\nsysdesign: Takeoff Runway Too Short, Increase Thrust')
    else:
        print('\nsysdesign: Takeoff Runway OK')
    if s > LandingRunwayLength:
        print('\nsysdesign: Landing Runway Too Short')
    else:
        print('\nsysdesign: Landing Runway OK')
    #prints advisement
    print('\nsysdesign: Rerun Designer To Ensure Accuracy')



# END OF ITERATIVE DESIGN SUBROUTINE


# START OF FLIGHT CALCULATIONS



elif RunSysDesign == False:


    # BEGINNING OF TAKE OFF CALCULATIONS
    
    
    
    #max lift coefficient assuming use of flaps
    CLmaxtakeoff = TakeoffLiftCoefficientWithFlaps*math.cos(math.radians(SweepAngle))
    #stall velocity at take off altitude
    Vstall = math.sqrt((2*(Weight/S))/(GetDensity(TakeoffAltitude)*CLmaxtakeoff))
    #calculates take off velocity
    Vtakeoff = 1.2*Vstall
    #takeoff acceleration
    atakeoff = 32.17*((EngineThrust*Engines)/Weight)
    #take off distance
    dtakeoff = pow(Vtakeoff,2)/(2*atakeoff)
    #take off time
    ttakeoff = Vtakeoff/atakeoff
    #take off thrust
    Thrusttakeoff = ((Weight/32.17)*math.sqrt(pow(Vtakeoff,2)/TakeoffRunwayLength))/Engines

    
    
    # END OF TAKE OFF CALCULATIONS
    
    
    # BEGINNING OF CLIMB CALCULATIONS
    
    
    
    #creates values for loop
    alt = TakeoffAltitude
    tclimb = 0
    Wclimb = 0
    newWeight = Weight
    #loop that calculates flight from take off location to cruise altitude
    while alt <= CruiseAltitude:
        #finds thrust at altitudes
        Talt1 = (GetDensity(alt)/GetDensity(0))*(EngineThrust*Engines)
        Talt2 = (GetDensity(alt+ClimbStep)/GetDensity(0))*(EngineThrust*Engines)
        #values that makes later calculations easier
        gamma1 = 1 + math.sqrt(1+((12*K*cd0)/pow((Talt1/newWeight),2)))
        gamma2 = 1 + math.sqrt(1+((12*K*cd0)/pow((Talt2/newWeight),2)))
        #fastest climb velocity
        Vfastestclimb1 = math.sqrt((((Talt1)/S)*gamma1)/(3*GetDensity(TakeoffAltitude)*cd0))
        Vfastestclimb2 = math.sqrt((((Talt2)/S)*gamma2)/(3*GetDensity(TakeoffAltitude)*cd0))
        #sine of the flight path angle
        singammamax1 = ((Talt1)/newWeight)*(1-(gamma1/6))-((6*K*cd0)/(gamma1*((Talt1)/newWeight)))
        singammamax2 = ((Talt2)/newWeight)*(1-(gamma2/6))-((6*K*cd0)/(gamma2*((Talt2)/newWeight)))
        #vertical velocity
        RoCmax1 = Vfastestclimb1*singammamax1
        RoCmax2 = Vfastestclimb2*singammamax2
        RoCavg = (RoCmax1+RoCmax2)/2
        #average thrust
        Tavg = (Talt1+Talt2)/2
        #time to climb
        dt = ClimbStep/RoCavg
        tclimb+=dt
        #fuel burned during climb
        dW = (SpecificFuelConsumption*(1/3600))*(Tavg*dt)
        #weight at end of step
        newWeight-=dW
        #weight loss during climb
        Wclimb+=dW
        #increments altitude
        alt+=ClimbStep
    #aircraft weight at start of cruise
    Wstartcruise = Weight - Wclimb
    
    #equations for fastest climb at sea level
    Talt = (GetDensity(TakeoffAltitude)/GetDensity(0))*(EngineThrust*Engines)
    gamma = 1 + math.sqrt(1+((12*K*cd0)/pow((Talt/Weight),2)))
    Vfastestclimb = math.sqrt((((Talt)/S)*gamma)/(3*GetDensity(TakeoffAltitude)*cd0))
    singammamaxfc = ((Talt)/Weight)*(1-(gamma/6))-((6*K*cd0)/(gamma*(Talt/Weight)))
    gammafc = math.degrees(math.asin(singammamaxfc))
    RoCfc = Vfastestclimb*singammamaxfc*60
    #equations for steepest climp at sea level
    Vsteepestclimb = math.sqrt(((2*(Weight/S))/GetDensity(TakeoffAltitude))*math.sqrt(K/cd0))
    singammamaxsc = (Talt/Weight)-(2*math.sqrt(K*cd0))
    gammasc = math.degrees(math.asin(singammamaxsc))
    RoCsc = Vsteepestclimb*singammamaxsc*60
    
    
    
    # END OF CLIMB CALCULATIONS
    
    
    # BEGINNING OF STEADY LEVEL FLIGHT CALCULATIONS
    
    
    
    #time of flight assuming constant altitude and constant velocity
    t = (FlightRange*5280)/CruiseAirspeed
    #estimated weight loss
    dW = SpecificFuelConsumption*(1/3600)*(EngineThrust*Engines*(GetDensity(CruiseAltitude)/GetDensity(0)))*t
    Wendcruise = Wstartcruise-dW
    
    #required thrust at altitude cruise
    Treqcruisestart = (0.5*GetDensity(CruiseAltitude)*pow(CruiseAirspeed,2)*cd0*S)+((K*pow(Wstartcruise/S,2)*S)/(0.5*GetDensity(CruiseAltitude)*pow(CruiseAirspeed,2)))
    Treqcruiseend = (0.5*GetDensity(CruiseAltitude)*pow(CruiseAirspeed,2)*cd0*S)+((K*pow(Wendcruise/S,2)*S)/(0.5*GetDensity(CruiseAltitude)*pow(CruiseAirspeed,2)))
    #average thrust during cruise
    Treqavg = (Treqcruisestart+Treqcruiseend)/2
    #average thrust per engine
    Treqavgeng = Treqavg/Engines
    #average weight during cruise
    Wavg = (Wstartcruise + Wendcruise)/2
    #average efficiency during flight
    Eavg = Wavg/Treqavg
    #best range velocity
    Vbr = math.sqrt(((2*(Wstartcruise/S))/GetDensity(CruiseAltitude))*math.sqrt((3*K)/cd0))
    #finds weight of wing
    Wwing = (0.0051*pow(Weight*3.84,0.557)*pow(S,0.649)*pow(AR,0.5)*pow(1+(TipChord/RootChord),0.1)*pow(S*0.1,0.1))/(pow(MaxThickness,0.4)*pow(math.cos(math.radians(SweepAngle)),1))
    
    
    
    # END OF STEADY LEVEL FLIGHT CALCULATIONS
    
    
    # BEGINNING OF LEVEL TURNING CALCULATIONS
    
    
    
    #sets maximum load factor
    nmax = 3.84
    #velocity for fastest turn
    Vfastestturn = math.sqrt(((2*(Wavg/S))/GetDensity(CruiseAltitude))*math.sqrt(K/cd0))
    #finds load factor for fastest turn
    Talt = EngineThrust*Engines*(GetDensity(CruiseAltitude)/GetDensity(0))
    n1 = math.sqrt(((Talt/Wavg)/math.sqrt(cd0*K))-1)
    if n1 < nmax:
        #for load factors under max allowable
        #load factor for fastest turn
        nft = n1
        #thrust for fastest turn
        Tft = Talt
    elif n1 >= nmax:
        #for load factors over max allowable
        #load factor for fastest turn
        nft = nmax
        #thrust for fastest turn
        Tft = Wavg*(pow(nft,2)+1)*math.sqrt(K*cd0)
    
    #yaw rate for fastest turn
    yawrateft = (32.17*math.sqrt(pow(nft,2)-1))/Vfastestturn
    #bankangle for fastest turn
    bankangleft = math.degrees(math.acos(1/nft))
    #fastest turn radius
    turnradiusft = pow(Vfastestturn,2)/(32.2*math.sqrt(pow(nft,2)-1))
    
    #load factor for minimum turn radius
    n2 = math.sqrt(2-((4*K*cd0)/pow(Talt/Wavg,2)))
    #two different stall airspeeds
    Vstall1 = math.sqrt(((2*(Wavg/S))/(GetDensity(CruiseAltitude)*CLmax))*n2)
    Vstall2 = math.sqrt((2*(Wavg/S)*nmax)/(GetDensity(CruiseAltitude)*CLmax))
    #airspeeds to compare to stall airspeed
    Vcomp1 = math.sqrt((4*K*(Wavg/S))/(GetDensity(CruiseAltitude)*(Talt/Wavg)))
    Vcomp2 = math.sqrt(((Wavg/S)*abs(2-pow(nmax,2)))/(GetDensity(CruiseAltitude)*cd0))
    if n2 < nmax:
        #for load factors under max allowable
        if Vcomp1 > Vstall1:
            #tightest turn thrust
            Ttt = Talt
            #tightest turn airspeed
            Vtightestturn = Vcomp1
            #tightest turn load factor
            ntt = n2
        elif Vcomp1 <= Vstall1:
            #tightest turn thrust
            Ttt = Wavg*math.sqrt(2*K*(cd0+(K*pow(CLmax,2))))
            #tightest turn load factor
            ntt = Clmax*math.sqrt((2*K)/(2*K*(cd0+(K*pow(CLmax,2)))))
            #tightest turn airspeed
            Vtightestturn = Vstall1
    elif n2 >= nmax:
        #for load factors over max allowable
        if Vcomp2 > Vstall2:
            #tightest turn load factor
            ntt = nmax
            #tightest turn thrust
            Ttt = (4*K*cd0*Wavg)/abs(2-pow(nmax,2))
            #tightest turn airspeed
            Vtightestturn = Vcomp2
        elif Vcomp2 <= Vstall2:
            #tightest turn load factor
            ntt = nmax
            #tightest turn airspeed
            Vtightestturn = Vstall2
            #tightest turn thrust
            Ttt = (2*K*Wavg*CLmax)/nmax
            
    #yaw rate for tightest turn
    yawratett = (32.17*math.sqrt(pow(ntt,2)-1))/Vtightestturn
    #bankangle for tightest turn
    bankanglett = math.degrees(math.acos(1/ntt))
    #tightest turn radius
    turnradiustt = pow(Vtightestturn,2)/(32.2*math.sqrt(pow(ntt,2)-1))    
    
    
    
    # END OF LEVEL TURNING FLIGHT CALCULATIONS
    
    
    # BEGINNING OF LANDING CALCULATIONS
    
    
    
    #finds max wing lift coefficient at landing
    CLmaxland = LandingLiftCoefficentWithFlaps*math.cos(math.radians(SweepAngle))
    #finds stall velocity at landing location
    Vstall = math.sqrt((2*(Wendcruise/S))/(GetDensity(LandingAltitude)*CLmaxland))
    #touchdown velocity
    Vtouchdown = 1.15*Vstall
    #distance traveled while rolling
    sfreeroll = Vtouchdown*PilotReactionTime
    #some things that make later calculations easier
    JT = (ReverseThrust*Engines)/Wendcruise
    JA = (GetDensity(LandingAltitude)*cd0)/(2*(Wendcruise/S))
    #finds distance required to stop
    sdecelerate = (1/(2*32.2*JA))*math.log(1+((JA/JT)*pow(Vtouchdown,2)),2.71828)
    s = sfreeroll + sdecelerate
    
    
    
    # END OF LANDING CALCULATIONS
    
    
    # BEGINNING OF UNPOWERED FLIGHT CALCULATIONS
    
    
    
    #best unpowered velocity at altitude and sea level
    Vbestunpowered = math.sqrt(((2*(Wavg/S))/(GetDensity(UnpoweredTestAltitude )))*math.sqrt(K/cd0))
    VbestunpoweredSL = math.sqrt(((2*(Wavg/S))/(GetDensity(0)))*math.sqrt(K/cd0))
    #max aerodynamic efficiency
    EM = (1/(2*math.sqrt(K*cd0)))
    #best unpowered range
    xunpowered = EM*(UnpoweredTestAltitude )
    #best unpowered flight time
    tunpowered = ((2*33333*EM)/VbestunpoweredSL)*(1-pow(2.71828,-UnpoweredTestAltitude /(2*33333)))
    #flight path angle
    gammaunpowered = math.degrees(math.atan(1/EM))
    
    
    
    # END OF UNPOWERED FLIGHT CALCULATIONS    
    
    
    # START OF ADVISEMENT CALCULATIONS
    
    
    
    #total fuel loss during flight
    Wlosstotal = Wclimb + dW
    #reserve fuel required
    Wreserve = SpecificFuelConsumption*(EngineThrust*Engines*(GetDensity(15000)/GetDensity(0)))*0.75
    #compares flight fuel loss to fuel available
    Weightdiff = WeightOfFuel - Wlosstotal
    #compares required fuel reserve to provided
    Weightreservediff = WeightOfFuelReserve - Wreserve
    
    #speed of sound at cruise altitude
    SoSadvice = math.sqrt(1.4*1716*GetTemperature(CruiseAltitude))
    #stall airspeed at cruise altitude
    Vstalladvice = math.sqrt((2*(Wavg/S))/(GetDensity(CruiseAltitude)*CLmax))
    #coffin airspeed at cruise altitude
    Vstaradvice = math.sqrt((2*(Wavg/S)*8)/(GetDensity(CruiseAltitude)*CLmax))
    #mach number for stall airspeed
    mstalladvice = Vstalladvice/SoSadvice
    #mach number for coffin airspeed
    mstaradvice = Vstaradvice/SoSadvice
    #mach number for cruise airspeed
    machadvice = CruiseAirspeed/SoSadvice
    
    #weight of all parts of aircraft
    We1 = WeightOfFuselage
    Wf = WeightOfFuel+WeightOfFuelReserve
    Wp = WeightOfCargo
    Wc = WeightOfPerson*People
    Wn = WeightOfEngines*Engines
    #total aircraft weight
    W0 = We1+Wf+Wp+Wc+Wn
    #estimate of empty weight using total weight
    We = EmptyWeightRatio*W0
    #compares estimate of empty weight to actual value
    Wediff = WeightOfFuselage-We
    
    #required thrust at altitude cruise
    Treqcruisestart = (0.5*GetDensity(CruiseAltitude)*pow(CruiseAirspeed,2)*cd0*S)+((K*pow(Wstartcruise/S,2)*S)/(0.5*GetDensity(CruiseAltitude)*pow(CruiseAirspeed,2)))
    Treqcruiseend = (0.5*GetDensity(CruiseAltitude)*pow(CruiseAirspeed,2)*cd0*S)+((K*pow(Wendcruise/S,2)*S)/(0.5*GetDensity(CruiseAltitude)*pow(CruiseAirspeed,2)))
    #average thrust during cruise
    Treqavg = (Treqcruisestart+Treqcruiseend)/2
    #average thrust per engine
    Treqavgeng = Treqavg/Engines
    #compares Thrust required to provided
    Treqdiff = Treqavgeng - EngineThrust
    
    
    
    # START OF OUTPUTS
    
    
    
    # Full length outputs
    
    
    if AbridgedOutputs == False:
        
        # Output for File Parsing and Aerodynamics
        
        
        #displays graph of airfoil lift curve
        print('sysprint: Airfoil Lift Curve')
        plt.subplot(121)
        plt.plot(apts,clpts)
        plt.title('C_l vs. Alpha (Airfoil)')
        plt.xlabel('Alpha (deg)')
        plt.ylabel('Coefficient of Lift')
        #displays graph of wing lift curve
        print('sysprint: Wing Lift Curve')
        plt.subplot(122)
        plt.plot(range(10),CLpts)
        plt.title('C_L vs Alpha (Wing)')
        plt.xlabel('Alpha (deg)')
        plt.ylabel('C_L')
        
        #prints section header
        print('\nAIRFOIL INFORMATION')
        #prints maximum lift coefficient
        print('Maximum Lift Coefficient:\n',round(CLmax,2))
        #prints AoA where lift is zero
        print('Zero Lift AoA:\n',round(a0,2),'deg') 
        #prints slope of airfoil lift curve
        print('Slope Of Airfoil Lift Curve:\n',round(a_0,3),'1/deg')
        #prints slope of wing lift curve
        print('Slope of Wing Lift Curve:\n',round(a,3),'1/deg')
        
        #prints section header
        print('\nWING INFORMATION')
        #prints wing area
        print('Wing Planform Area:\n',round(S,2),'ft^2')
        #prints average chord length
        print('Average Chord Length:\n',round(chord,2),'ft')
        #prints taper ratio
        print('Taper Ratio:\n',round(TipChord/RootChord,2))
        #prints wing zero lift coefficient
        print('Wing Zero Lift Drag:\n',round(cd0wing,3))
        #prints oswald efficiency factor
        print('Oswald Efficiency Factor:\n',round(OEF,2))
        #prints wing mount angle
        print('Wing Mounting Angle:\n',round(mountang,2),'deg')
        #prints lift coeficient for max aerodynamic efficiency
        print('Lift Coefficient for Max Efficiency:\n',round(CLEM,2))
        #prints angle of attack for max aerodynamic efficiency
        print('AoA for Max Efficiency:\n',round(aEM,2),'deg')
        #prints wing weight
        print('Total Wing Weight:\n',round(Wwing,2),'lb')
        
        #prints section header
        print('\nFUSELAGE INFORMATION')
        #prints wetted surface area of fuselage
        print('Wetted Fuselage Surface:\n',round(Swet,2),'ft^2')
        #prints fuselage zerolift drag coefficient
        print('Fuselage Zero Lift Drag:\n',round(cd0fuselage,3))
        
        #prints section header
        print('\nAIRCRAFT INFORMATION')
        #prints drag polar for aircraft
        print('Drag Polar:')
        print(' C_D = {} + {}*C_L^2'.format(round(cd0,4),round(K,4)))
        #prints aircrafts total weight
        print('Total Aircraft Weight:\n',round(Weight,2),'lb')
        #prints maximum efficiency
        print('Maximum Efficiency:\n',round(Emax,2))
        #prints zero lift drag for engine nacelle
        print('Nacelle Zero Lift Drag:\n',round(cd0nacelle,3))
        
        
        # Take off Output:
        
        
        #prints section header
        print('\nTAKEOFF INFORMATION')
        #prints takeoff wing loading
        print('Takeoff Wing Loading:\n',round(Weight/S,2),'psf')
        #prints take off velocity
        print('Takeoff Velocity:\n',round(Vtakeoff,2),'fps')
        #prints take off distance
        print('Takeoff Distance:\n',round(dtakeoff,2),'ft')
        #prints time for take off
        print('Takeoff Time:\n',round(ttakeoff,2),'sec')
        #prints take off thrust
        print('Thrust For Takeoff Per Engine:\n',round(Thrusttakeoff,2),'lb')
        #prints best range airspeed
        print('Best Range Airspeed:\n',round(Vbr,2),'fps')
        
        
        # Climb Outputs:
        
        
        #prints section header
        print('\nCLIMB INFORMATION')
        #angle for steepest climb
        print('AoA For Steepest Climb:\n',round(gammasc,2),'deg')
        #airspeed for steepest climb
        print('Steepest Climb Airspeed:\n',round(Vsteepestclimb,2),'fps')
        #steepest rate of climb
        print('Steepest Rate Of Climb:\n',round(RoCsc,2),'fpm')
        #angle for fastest climb
        print('AoA For Fastest Climb:\n',round(gammafc,2),'deg')
        #airspeed for fastest climb
        print('Fastest Climb Airspeed:\n',round(Vfastestclimb,2),'fps')
        #fastest rate of climb
        print('Fastest Rate Of Climb:\n',round(RoCfc,2),'fpm')
        #prints time taken to climb
        print('Time to Climb:\n',round(tclimb/60,2),'min')    
        
        
        # Output for steady level flight
        
        
        #prints section header
        print('\nCRUISE INFORMATION')
        #prints total flight time
        print('Flight Time:\n',round(t/3600,2),'hr')
        #prints cruise required thrust
        print('Thrust Required For Cruise Per Engine:\n',round(Treqavgeng,2),'lb')
        #prints efficiency during cruise
        print('Cruise Aerodynamic Efficiency:\n',round(Eavg,2))
        #prints weight loss during cruise
        print('Weight Of Fuel Burned:\n',round(dW,2),'lb')
        #prints weight at end of cruise
        print('Weight After Cruise:\n',round(Wendcruise,2),'lb')
        
        
        # Outputs for level turning flight
        
        
        #prints section header
        print('\nTURNING FLIGHT INFORMATION')
        #prints fastest turn bank angle
        print('Fastest Turn Bank Angle:\n',round(bankangleft,2),'deg')
        #prints fastest turn load factor
        print('Fastest Turn Load Factor:\n',round(nft,2))
        #prints thrust for fastest turn
        print('Fastest Turn Thrust Per Engine:\n',round(Tft/Engines,2),'lb')
        #prints airspeed for fastest turn
        print('Fastest Turn Airspeed:\n',round(Vfastestturn,2),'fps')
        #prints turn radius for fastest turn
        print('Fastest Turn Radius:\n',round(turnradiusft,2),'ft')
        #prints fastest turn yaw rate
        print('Fastest Turn Yaw Rate:\n',round(yawrateft,2),'rad/s')
        
        #prints tightest turn bank angle
        print('Tightest Turn Bank Angle:\n',round(bankanglett,2),'deg')
        #prints tightest turn load factor
        print('Tightest Turn Load Factor:\n',round(ntt,2))
        #prints thrust for tightest turn
        print('Tightest Turn Thrust Per Engine:\n',round(Ttt/Engines,2),'lb')
        #prints airspeed for tightest turn
        print('Tightest Turn Airspeed:\n',round(Vtightestturn,2),'fps')
        #prints turn radius for tightest turn
        print('Tightest Turn Radius:\n',round(turnradiustt,2),'ft')
        #prints tightest turn yaw rate
        print('Tightest Turn Yaw Rate:\n',round(yawratett,2),'rad/s')
        
        
        # Outputs for Landing
        
        
        #prints section header
        print('\nLANDING INFORMATION')
        #prints touchdown velocity
        print('Touchdown Velocity:\n',round(Vtouchdown,2),'fps')
        #prints JA
        print('JA:\n',round(JA,9))
        #prints JT
        print('JT:\n',round(JT,2))
        #prints free roll distance
        print('Free Roll Distance:\n',round(sfreeroll,2),'ft')
        #prints distance to decelerate
        print('Deceleration Distance:\n',round(sdecelerate,2),'ft')
        #prints distance required to stop
        print('Ground Roll Distance:\n',round(s,2),'ft')
        #prints landing wing loading
        print('Landing Wing Loading:\n',round(Wendcruise/S,2),'psf')
        #prints weight at landing
        print('Landing Weight:\n',round(Wendcruise,2),'lb')
        
        
        #Outputs for Unpowered Flight
        
        
        #prints section header
        print('\nUNPOWERED FLIGHT INFORMATION')
        #prints best range
        print('Best Unpowered Range:\n',round(xunpowered/5280,2),'mi')
        #prints velocity for best range
        print('Velocity for Best Unpowered Range:\n',round(Vbestunpowered,2),'fps')
        #prints endurance for best range
        print('Endurance for Best Unpowered Range:\n',round(tunpowered/60,2),'min')
        #prints flight angle for best range
        print('Flight Path Angle:\n',round(gammaunpowered,2),'deg')
        
        
        
    # Abridged Outputs
    
    
    elif AbridgedOutputs == True:
        
        #prints drag polar for aircraft
        print('Drag Polar:')
        print(' C_D = {} + {}*C_L^2'.format(round(cd0,4),round(K,4)))
        #prints aircrafts total weight
        print('Total Aircraft Weight:\n',round(Weight,2),'lb')
        #prints maximum efficiency
        print('Maximum Efficiency:\n',round(Emax,2))
        #prints take off velocity
        print('\nTakeoff Velocity:\n',round(Vtakeoff,2),'fps')
        #prints take off distance
        print('Takeoff Distance:\n',round(dtakeoff,2),'ft')
        #prints time for take off    
        print('Takeoff Time:\n',round(ttakeoff,2),'sec')
        #prints time to climb
        print('\nTime to Climb:\n',round(tclimb/60,2),'min')    
        #prints total flight time
        print('\nFlight Time:\n',round(t/3600,2),'hr')
        #prints cruise required thrust
        print('Thrust Required For Cruise Per Engine:\n',round(Treqavgeng,2),'lb')
        #prints efficiency during cruise
        print('Cruise Aerodynamic Efficiency:\n',round(Eavg,2))
        #prints touchdown velocity
        print('\nTouchdown Velocity:\n',round(Vtouchdown,2),'fps')
        #prints distance required to stop
        print('Ground Roll Distance:\n',round(s,2),'ft')        
        
        
        
    # Error Statement For Output Settings
    
    
    else:
        print('\nsyserror: Incorrect Output Settings, Try Again')
        exit()
    
    
    
    # Outputs for Advisement Calculations
    
    
    #advisor system for fuselage weight
    if abs(Wediff) > 10:
        print('\nsysadvice: Recheck Fuselage Weight')
    else:
        print('\nsysadvice: Fuselage Weight OK')
    
    
    #advisor system for fuel weight
    if abs(Weightdiff) > 10:
        print('sysadvice: Recheck Fuel')
    else:
        print('sysadvice: Fuel OK')
    
    
    #advisor system for reserve fuel weight
    if abs(Weightreservediff) > 10:
        print('sysadvice: Recheck Reserve Fuel')
    else:
        print('sysadvice: Reserve Fuel OK')
    
    
    #advisor system for cruise speed
    if machadvice >= mstalladvice and machadvice <= mstaradvice:
        #runs if given cruise speed is within +- 0.05 of required
        print('sysadvice: Cruise Speed OK')
    elif machadvice < mstalladvice:
        #runs if given cruise speed is at most 0.05 too slow
        print('sysadvice: Increase Cruise Speed to',round(SoSadvice*0.8,2),'fps')
    elif machadvice > mstaradvice:
        #runs if given cruise speed is at least 0.05 too fast
        print('sysadvice: Decrease Cruise Speed to',round(SoSadvice*0.8,2),'fps')
    else:
        print('sysadvice: Error in Cruise Speed Calculations')
    
    
    #advisor system for Engine Thrust
    if Treqdiff <-10:
        print('sysadvice: Thrust OK')
    elif Treqdiff >= -10 and Treqdiff <= 10:
        print('sysadvice: Thrust OK')
    elif Treqdiff > 10:
        print('sysadvice: Consider New Thrust',round(Treqavgeng,2),'lb')
    else:
        print('sysadvice: Error In Thrust Calculations')
    
    
    
    # END OF OUTPUTS
    
    
    #This must be after all calculations
    plt.show()



# END OF FLIGHT CALCULATIONS


# ERROR CONDITION FOR DESIGN SUBROUTINE MODES



#runs if a boolean is not entered for 'RunSysDesign'
else:
    print('syserror: Incorrect Input For Design Controls, Try Again')



# END OF ERROR CONDITION


# END OF PROGRAM