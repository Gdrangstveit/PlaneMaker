# START OF PROGRAM


#imports math library and plotting library
import math
import matplotlib.pyplot as plt


# PROGRAM INFORMATION


#PLANE MAKER 2
#CREATED BY GUNNAR DRANGSTVEIT
#v<2.3.1>

#Read the ReadMe file that came with this program for instructions and related information.
#Design Mode: Change 'RunSysDesign' to True if you wish to enter Design mode, and False if you do not


# USER INPUTS (ENTER ALL IN DESIGNATED IMPERIAL UNITS):
#(initial inputs are for the Boeing 747, change all and delete this line)

#fuselage dimensions
FuselageLength = 250 #ft
FuselageDiameter = 18 #ft

#wing dimensions
RootChord = 30 #ft
TipChord = 6 #ft
WingSpan = 211 #ft
SweepAngle = 37.5 #deg
LocationOfMaxThickness = 0.35 #unitless
MaxThickness = 0.113 #unitless

#file containing lift coifficients and AoA's
FileName = 'C:/Users/user1/Documents/Scripts/Airfoil.txt' #Change instances of '\' to '/'

#cruise phase informations
CruiseAltitude = 30000 #ft
CruiseAirspeed = 123 #fps
FlightRange = 6185 #mi

#aircraft weight information
WeightOfFuselage = 398800 #lb
WeightOfPerson = 170 #lb
WeightOfCargo = 110000 #lb
WeightOfFuel = 162580 #lb
WeightOfFuelReserve = 0 #lb
WeightOfEngines = 8825 #lb

#onboard quantities
People = 600 #people
Engines = 4 #unitless (number of engines)

#engine information
EngineThrust = 54000 #lb (per engine)
SpecificFuelConsumption = 0.605 #1/hr
ReverseThrust = 15000 #lb

#information of takeoff, climb, and landing locations
TakeoffAltitude = 0 #ft
TakeoffRunwayLength = 12345 #ft
ClimbStep = 1000 #ft (sets accuracy of climb sim)
UnpoweredTestAltitude = 30000 #ft
LandingAltitude = 0 #ft
LandingRunwayLength = 12345 #ft
PilotReactionTime = 3 #sec
LandingLiftCoefficentWithFlaps = 2.54 #unitless (Found in Read Me)
TakeoffLiftCoefficientWithFlaps = 1.90 #unitless (found in Read Me)

#iterative design mode controls
RunSysDesign = False #boolean (True or False)
WeightAccuracy = 0.001 #lb (lower than 0.001 may take long run times)
EmptyWeightRatio = 0.5 #unitless (Found in Read Me)


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
    k=k+1
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
#form factor for wing
FFwing = (1+(0.6/LocationOfMaxThickness)*(MaxThickness)+100*pow(MaxThickness,4))*(1.34*pow(mach,0.18)*pow(math.cos(math.radians(SweepAngle)),0.28))
#form factor of fuselage
f = FuselageLength/FuselageDiameter
FFfuselage = 0.9 + (5/pow(f,1.5))+(f/400)
#wing area
S = WingSpan*chord
#surface area of fuselage
SAcone = 3.14159*(FuselageDiameter/2)*((FuselageDiameter/2)+math.sqrt(pow((FuselageLength*0.1),2)+pow(FuselageDiameter/2,2)))
SAcylinder = 2*3.14159*(FuselageDiameter/2)*(FuselageLength*0.9)
Swet = SAcone + SAcylinder
#no lift drag
cd0 = (cfwing*FFwing*2*2)+(cffuselage*FFfuselage*1*(Swet/S))

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


# END OF FILE PARSING AND AERODYNAMIC CALCULATIONS

# START OF ITERATIVE DESIGN SUBROUTINE


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
        CLmaxland = LandingLiftCoefficientWithFlaps*math.cos(math.radians(SweepAngle))
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
    # Output for File Parsing and Aerodynamics
    
    #displays graph of airfoil lift curve
    print('sysprint: Airfoil Lift Curve')
    plt.subplot(121)
    plt.plot(apts,clpts)
    plt.title('C_l vs. Alpha (Airfoil)')
    plt.xlabel('Alpha (deg)')
    plt.ylabel('Coefficient of Lift')
    #displays graph of wing lift curve
    print('sysprint: Wing Lift Curve\n')
    plt.subplot(122)
    plt.plot(range(10),CLpts)
    plt.title('C_L vs Alpha (Wing)')
    plt.xlabel('Alpha (deg)')
    plt.ylabel('C_L')
    #prints drag polar for aircraft
    print('Drag Polar:')
    print(' C_D = {} + {}*C_L^2'.format(round(cd0,4),round(K,4)))
    #prints aircrafts total weight
    print('Total Aircraft Weight:\n',round(Weight,2),'lb')
    #prints maximum efficiency
    print('Maximum Efficiency:\n',round(Emax,2))

    
    
    # BEGINNING OF TAKE OFF CALCULATIONS
    
    

    #max lift coefficient assuming use of flaps
    CLmaxtakeoff = TakeoffLiftCoefficentWithFlaps*math.cos(math.radians(SweepAngle))
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
    
    # Take off Output:
    
    #prints take off velocity
    print('\nTake off Velocity:\n',round(Vtakeoff,2),'fps')
    #prints take off distance
    print('Take off Distance:\n',round(dtakeoff,2),'ft')
    #prints time for take off
    print('Take off Time:\n',round(ttakeoff,2),'sec')
    
    
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
    
    # Climb Outputs:
    
    #prints time taken to climb
    print('\nTime to Climb:\n',round(tclimb/60,2),'min')
    
    
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
    
    # Output for steady level flight
    
    #prints total flight time
    print('\nFlight Time:\n',round(t/3600,2),'hr')
    #prints cruise required thrust
    print('Thrust Required For Cruise Per Engine:\n',round(Treqavgeng,2),'lb')
    #prints efficiency during cruise
    print('Cruise Aerodynamic Efficiency:\n',round(Eavg,2))
    
    
    # END OF STEADY LEVEL FLIGHT CALCULATIONS
    
    # BEGINNING OF LEVEL TURNING CALCULATIONS
    
    
    #velocity of minimum turn radius
    Vminradius = math.sqrt((4*K*(Wavg/S))/(GetDensity(CruiseAltitude)*((EngineThrust*Engines*(GetDensity(CruiseAltitude)/GetDensity(0)))/Wavg)))
    #load factor for minimum turn radius
    nminradius2 = 2-((4*K*cd0)/pow((EngineThrust*Engines*(GetDensity(CruiseAltitude)/GetDensity(0)))/Wavg,2))
    #minimum turn radius
    minradius = pow(Vminradius,2)/(32.2*math.sqrt(nminradius2-1))
    
    # Outputs for level turning flight
    
    #prints minimum turn radius
    print('\nMinimum Turn Radius:\n',round(minradius,2),'ft')
    #prints velocity for minimum turn radius
    print('Velocity for Minimum Turn Radius:\n',round(Vminradius,2),'fps')
    
    
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
    s = sfreeroll + (1/(2*32.2*JA))*math.log(1+((JA/JT)*pow(Vtouchdown,2)),2.71828)
    
    # Outputs for Landing
    
    #prints touchdown velocity
    print('\nTouchdown Velocity:\n',round(Vtouchdown,2),'fps')
    #prints distance required to stop
    print('Ground Roll Distance:\n',round(s,2),'ft')
    
    
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
    
    #Outputs for Unpowered Flight
    
    #prints best range
    print('\nBest Unpowered Range:\n',round(xunpowered/5280,2),'mi')
    #prints velocity for best range
    print('Velocity for Best Unpowered Range:\n',round(Vbestunpowered,2),'fps')
    #prints endurance for best range
    print('Endurance for Best Unpowered Range:\n',round(tunpowered/60,2),'min')
    #prints flight angle for best range
    print('Flight Path Angle:\n',round(gammaunpowered,2),'deg')
    
    
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
    We = 0.5*W0
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
    
    #advisor system for Engine Thrust
    if Treqdiff <-10:
        print('sysadvice: Thrust OK')
    elif Treqdiff >= -10 and Treqdiff <= 10:
        print('sysadvice: Thrust OK')
    elif Treqdiff > 10:
        print('sysadvice: Consider New Thrust',round(Treqavgeng,2),'lb')
    
    # END OF ADVISEMENT CALCULATIONS
    
    
    #This must be after all calculations
    plt.show()


# END OF FLIGHT CALCULATIONS

# ERROR CONDITION FOR DESIGN SUBROUTINE MODES


#runs if a boolean is not entered for 'RunSysDesign'
else:
    print('syserror: Incorrect Input For Design Controls, Try Again')


# END OF ERROR CONDITION


# END OF PROGRAM