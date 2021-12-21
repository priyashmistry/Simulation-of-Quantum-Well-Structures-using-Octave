clc;
clear all;
#-------------------------------------------------------------------------------
#                     Simulation of Quantum Confinement
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                                 Instructions
#-------------------------------------------------------------------------------
disp("-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--.-.-.-.-.-.-.-.-.-.-.-.")
disp("           WELCOME TO THE SIMULATION OF QUANTUM CONFINEMENT              ")
disp("-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--.-.-.-.-.-.-.-.-.-.-.-.")
disp("Here are some examples of Quantum Well Structure if you want to try ...  ")
disp(" ")
disp("1) GaAs & AlGaAs")
disp("  # Active Material : GaAs ")
disp("  # Cladding Material : AlGaAs ")
disp("  # Bandgap of Active Region : 1.42 eV ")
disp("  # Bandgap of Cladding Region : 1.80 eV ")
disp("  # Effective Mass of Electron : 0.063 ")
disp("  # Effective Mass of Light Hole : 0.082 ")
disp("  # Effective Mass of Heavy Hole : 0.54 ")
disp("  # Thermal Velocity of Electron : 4.4 * 10^5 m/s ")
disp(" ")
disp("2) GaInAs & InP ")
disp("  # Active Material : GaInAs ")
disp("  # Cladding Material :  InP ")
disp("  # Bandgap of Active Region : 0.75 eV ")
disp("  # Bandgap of Cladding Region : 1.35 eV ")
disp("  # Effective Mass of Electron : 0.041 ")
disp("  # Effective Mass of Light Hole : 0.051 ")
disp("  # Effective Mass of Heavy Hole : 0.47 ")
disp("  # Thermal Velocity of Electron : 5.5 * 10^5 m/s ")
disp(" ")
disp("3) InGaAs & GaAs ")
disp("  # Active Material : InGaAs ")
disp("  # Cladding Material : GaAs ")
disp("  # Bandgap of Active Region : 0.75 eV ")
disp("  # Bandgap of Cladding Region : 1.42 eV ")
disp("  # Effective Mass of Electron : 0.041 ")
disp("  # Effective Mass of Light Hole : 0.082 ")
disp("  # Effective Mass of Heavy Hole : 0.54 ")
disp("  # Thermal Velocity of Electron : 5.5 * 10^5 m/s ")
disp(" ")
disp("-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.--.-.-.-.-.-.-.-.-.-.-.-.")
disp(" ")
disp(" ")
disp(".........................................................................")
#-------------------------------------------------------------------------------
#                Want to try from above or your own material
#-------------------------------------------------------------------------------

option1 = input("Do you want to try for one of the above Materials ?? (1 = Yes, 0 = No) : ");

# Want to try from above material
if option1==1
  yes = input("Enter the number 1,2 or 3 : ");
  
  # 1st material
  if yes==1 
    disp("You chose the material : GaAs/AlGaAs ")
    disp(".........................................................................")
    Eg = 1.42;
    Eg2 = 1.80;
    Me = 0.063;
    Ml = 0.082;
    Mh = 0.54;  
    V = 4.4*10^5;
    
  # 2nd material
  elseif yes==2 
    disp("You chose the material : GaInAs/InP ")
    disp(".........................................................................")
    Eg = 0.75;
    Eg2 = 1.35;
    Me = 0.041;
    Ml = 0.051;
    Mh = 0.47;
    V = 5.5*10^5;
    
  # 3rd Material  
  elseif yes==3 
    disp("You chose the material : GaAs/InGaAs ")
    disp(".........................................................................")
    Eg = 0.75;
    Eg2 = 1.42;
    Me = 0.041;
    Ml = 0.082;
    Mh = 0.54;
    V = 5.5*10^5;
  
  # Incorrect Value
  else 
    disp("Please Enter Correct Value")
    return
  endif
  Lo = input("Please enter width of Quantum Well(in nm) : ");

# Want to try your own material
elseif option1==0
disp(".........................................................................")
disp("                  Input the Values for your own material                 ")
disp(".........................................................................")
Eg = input("Please enter bandgap of active region in eV : ");
Eg2 = input ("Bandgap of cladding material in eV : ");
Me = input("Effective mass of electron in 'x*M0' form : ");
Ml = input("Effective mass of light holes in 'x*M0' form : "); 
Mh = input("Effective mass of heavy holes in 'x*M0' form : "); 
V = input("Thermal velocity of electrons : ");
Lo = input("Please enter width of Quantum Well(in nm) : ");

# Incorrect Value
else 
disp("Please Enter Correct Value")
return

endif 
# Average effective mass of holes
Mv = ((Ml)^1.5 + (Mh)^1.5)^(2/3); 

L = (0.6069/(sqrt(Eg2 - Eg)))*sqrt(1/Me + 1/Ml);  
disp(".........................................................................")
display("For EK Diagram")
disp(" ")
Efc = input("Please Enter Fermi Conduction Bandgap of Material in eV : ");
Efv = input("Please Enter Fermi Valence Bandgap of Material in eV : ");
disp(".........................................................................")


#-------------------------------------------------------------------------------
#                        EK Diagram for Active Region 
#-------------------------------------------------------------------------------
Evm = []; # Valence Bandgap Energy Array
Ecm = []; # Conduction Bandgap Energy Array
Ev1m = []; # Light Hole Valence Bandgap Energy Array
Ev2m = []; # Heavy Hole Valence Bandgap Energy Array
Km = []; # Wave-vector Array
h = 1*10^(-34); # Reduced Plank Constant
n = -10; 
for i=1:1:21
  n++;
  K = n*pi/(Lo*10^(-9)); # Wave Vector

  Ev1 = Efv - ((h^2)*(K^2))/(2*Ml*9.31*10^(-31)*1.6*10^(-19)); # Energy of light holes Valence Band  
  Ev2 = Efv - ((h^2)*(K^2))/(2*Mh*9.31*10^(-31)*1.6*10^(-19)); # Energy of heavy holes Valence Band 
  Ev = Efv - ((h^2)*(K^2))/(2*Mv*9.31*10^(-31)*1.6*10^(-19)); # Energy of Valence Band
  Ec = Efc + ((h^2)*(K^2))/(2*Me*9.31*10^(-31)*1.6*10^(-19)); # Energy of Conduction Band
  Evm(i) = Ev; # Store the data in array
  Ecm(i) = Ec; # Store the data in array
  Ev1m(i) = Ev1; # Store the data in array
  Ev2m(i) = Ev2; # Store the data in array
  Km(i) = K; # Store the data in array
  i++;
end
#-------------------------------------------------------------------------------
#                            Figure (1) : EK Diagram 
#-------------------------------------------------------------------------------
figure(1)
  
  plot(Km,Evm,"linewidth",3,'m');
  hold on;
  plot(Km,Ecm,"linewidth",3,'k');
  hold on
  plot(Km,Ev1m,"linewidth",3,'--y');
  hold on;
  plot(Km,Ev2m,"linewidth",3,'--g');
  hold on
  
legend("Valence Band (Average Effective Mass of Holes)","Conduction Band","Valence Band (Light Holes)","Valence Band (Heavy Holes)",
       "location","south")
  
  str1=sprintf("Bandgap E_g = %d eV \n",Efc-Efv);
  text(0,(Efc+Efv)/2,{str1},"horizontalalignment","center",'FontSize',15,'Color',[0.2,0.4,1],'FontWeight','bold') 
  str2=sprintf(" ");
  str3=sprintf("E-K Diagram for Active Material");
  set(gca,"linewidth",2,"fontsize",12,'FontWeight','bold')
  title({str2,str3,str2},'FontSize',20,'Color',[0,0,0],'FontWeight','bold')
  
  xlabel("Wave Vector k (m)",'FontSize',13,'Color','k','FontWeight','bold')
  ylabel("Energy E (eV)",'FontSize',13,'Color','k','FontWeight','bold')
  
  axis("square")
 

 
#-------------------------------------------------------------------------------
#              Check Whether Quantum Confinement is Possible or not 
#-------------------------------------------------------------------------------
Lamda_D = 6.63*10^(-25)/(Me*V*9.31*10^(-31)); 
# Condition : 1 when quantum confinement will not observe 
if ( Lo > Lamda_D)
  Lamda_E = 1.23*1000/Eg ; # Maximum Wavelength that can be emitted
  printf("Emission Wavelenth = %d \n",Lamda_E);
  printf("\n Effects of Quantum Confinement not observed. \n");
  printf("\n Width of Quantum Well should be equal to or less than De-Broglie Wavelength of Electron \n");  
  printf("\n De-Broglie Wavelength of electron is = %d \n",Lamda_D);
  return;

# Condition : 2 when quantum confinement will not observe 
elseif ( Lo < L)
  printf("For given material at %d nm (Quantum Well Width), Length of well exceeds minimum practical length \n",Lo)
  printf("Please enter length of well that exceeds > %d nm\n",L)
  return;
  
# Condition : 3 when quantum confinement will observe 
else 
  Lamda_E = 1.23*1000/(Eg + (0.369/Lo^2)*(1/Me + 1/Mv)); # Emission Wavelenght
  display("\n Here are your results");
  printf("\n # Emission Wavelenth = %d nm\n",Lamda_E);
endif

Lamda_max = 1.23*1000/Eg; # Maximum Emission Wavelength
Lamda_Min = 1.23*1000/Eg2; # Minimum Emission Wavelength

printf("\n # Practical Minimum Emission Wavelenth = %d nm\n",Lamda_Min);
printf("\n # Practical Maximum Emission Wavelenth = %d nm\n",Lamda_max);
printf("\n # Minimum optimum well width = %d nm\n",L);
printf("\n # Maximum well width upto which Quantum Confinment can be observed = %d nm\n",Lamda_D);
disp(".........................................................................")

#-------------------------------------------------------------------------------
#      Calculation of Wavelength as a function of Decreasing Well Width 
#-------------------------------------------------------------------------------
Lamda_E1_matrix = []; # Emission Wavelength Array 
Lm=[]; # Width of Well Array
counter=1; # Counter

for i = L:0.08:Lo
  Lm(counter)=i;
  Lamda_E1 = 1.23*1000/(Eg + (0.369/i^2)*(1/Me + 1/Mv)); # Emission Wavelength with Decreasing Width
  
  # If wavelength will be higher than the maximum it will become constant to it max value 
  if Lamda_E1>Lamda_max
    Lamda_E1=Lamda_max;
  
  # If wavelength will be lower than the minimum it will become constant to it min value   
  elseif Lamda_E1<Lamda_Min
    Lamda_E1=Lamda_Min;
    
  endif
  
  Lamda_E1_matrix(counter) = Lamda_E1; # Store the wavelength value in array 
  counter+=1;
end
#-------------------------------------------------------------------------------
#             Figure (2) :  Wavelength vs Width of Well Graph
#-------------------------------------------------------------------------------
figure(2)

plot(Lm,Lamda_E1_matrix,'r',"linewidth",3);

ylabel("Emission Wavelength (nm)",'FontSize',13,'Color','k','FontWeight','bold');
xlabel("Length of Well (nm)",'FontSize',13,'Color','k','FontWeight','bold');

str2=sprintf(" ");
str3=sprintf("Wavelength vs Width of Well");
set(gca,"linewidth",2,"fontsize",12,'FontWeight','bold')
title({str2,str3,str2},'FontSize',20,'Color',[0,0,0],'FontWeight','bold')

axis("square")

 
 
#-------------------------------------------------------------------------------
#                      Calculations to create Animation
#-------------------------------------------------------------------------------
v=20; # Initial Velocity of particle

theta_d = 23; # Initial angle of particle in xy direction
theta=theta_d*pi/180; # Convert theta in radian

phi_d = 23; # Initial angle of particle in z direction
phi=phi_d*pi/180; # Convert theta in radian

vx=v*cos(theta)*sin(phi); # Calculate x component of velocity
vy=v*sin(theta)*sin(phi); # Calculate y component of velocity
vz=v*cos(phi); # Calculate z component of velocity

x=2; # Initial x coordinate of particle 
y=3; # Initial y coordinate of particle 
z=1; # Initial z coordinate of particle 

a = 30; # Acceleration



#-------------------------------------------------------------------------------
#Figure (3) :  Animation How Electron moves in Quantum Well of Decreasing Width 
#-------------------------------------------------------------------------------
figure(3)

# Calculations for moving particle
for i=1:counter-1
  
x = x + vx*0.1;
y = y + vy*0.1;
z = z + vz*0.1;

vx = vx + a*0.1;
vy = vy + a*0.1;
vz = vz + a*0.1;

if x<0
  x=0;
  vx=-vx;
  
elseif x>=Lm(counter-i)
  x = Lm(counter-i);
  vx = -vx;
  
endif


if y<0
  y = 0;
  vy = -vy;
  
elseif y>=Lo
  y = Lo;
  vy = -vy;
  
endif 


if z<0
  z = 0;
  vz = -vz;
  
elseif z>=Lo
  z = Lo;
  vz = -vz;
  
endif 

# Creating 3D Box of decrasing width 
line_0toLo = linspace(0,Lo,counter-1); 
line_xyz0 = zeros(1,counter-1);
line_fixLo = Lo*ones(1,counter-1);

plot3(line_xyz0,line_xyz0,line_0toLo,'k','linewidth',5);
hold on 
plot3(line_xyz0,line_fixLo,line_0toLo,'k','linewidth',5);
hold on 
plot3(line_xyz0,line_0toLo,line_fixLo,'k','linewidth',5);
hold on 

q=quiver3(0,0,Lo,Lm(counter-i),0,0,'k','linewidth',5,'ShowArrowHead','off');
r=quiver3(0,Lo,Lo,Lm(counter-i),0,0,'k','linewidth',5,'ShowArrowHead','off');

s=quiver3(Lm(counter-i),0,0,0,0,Lo,'k','linewidth',5,'ShowArrowHead','off');
t=quiver3(Lm(counter-i),Lo,0,0,0,Lo,'k','linewidth',5,'ShowArrowHead','off');
u=quiver3(Lm(counter-i),0,Lo,0,Lo,0,'k','linewidth',5,'ShowArrowHead','off');

b=rectangle('Position',[0 0 Lm(counter-i) Lo],'linewidth',5);
hold on 
c=plot3(x,y,z,'.r','markersize',50);



str2=sprintf(" ");
str3=sprintf("Electron moving in Quantum Well of Decreasing Width");
set(gca,"linewidth",2,"fontsize",12,'FontWeight','bold')
title({str2,str3,str2},'FontSize',20,'Color',[0,0,0],'FontWeight','bold')

xlabel("Width of Quantum Well",'FontSize',13,'Color','k','FontWeight','bold')
axis([0 Lo 0 Lo 0 Lo],"square")

pause(0.01) 

if i<(counter-1)
  delete(c)
  delete(b)
  delete(q)
  delete(r)
  delete(s)
  delete(t)
  delete(u)
endif
 
end 

#-------------------------------------------------------------------------------
#                                 END OF PROGRAM 
#-------------------------------------------------------------------------------
