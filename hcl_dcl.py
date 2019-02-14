''' Name = Marty Krakora
Date = 13 February 2019
Title = HCl DCl analysis '''

from math import *

# After fitting data points to curve, fill in values here to get the constants of B_0, B_1, v_1, etc.
#Note here we are detecting the overtone transition so m = 2
#for each isotope, inserted the correct curve fit as well as the molar mass
x = 35 #molar mass of Cl
y = 1 #molar mass of H
v_m= "5674 + 20.3x -.712x^2"
delta_v = "2m(B_1 - B_0) + 2B_1"
a = -.712
b = 20.3
c_0 = 5674

#constants
c_speed = 3*(10**8) # m/s
h_bar = (((3.16152649*(10**(-26)))*(5.034*(10**22))) / (100)) #should be unitless
def mu(x, y):
    return ((x/1000)*(y/1000))/((x/1000)+(y/1000))

#calculations
v_0= 2885.9
v_1 = c_0
alpha_e = -(a) #units are wavenumbers
B_0 = .5*(b+(2*alpha_e)) #units are wavenumbers
B_1 = .5*b #units are wavenumbers
r_0 = (h_bar / (4*pi*mu(x,y)*B_0*c_speed))**(.5)
r_1 = (h_bar / (4*pi*mu(x,y)*B_1*c_speed))**(.5)
r_e = 1.267460*(10**(-10))
B_e = (h_bar / (4*pi*mu(x,y)*(r_e**2)*c_speed))
x_e = (2*v_0-v_1)/(2*v_0 - 2*v_1) #unitless
w_e = v_0 / (1-2*x_e) #units of wavenumbers
k = ((2*pi*c_speed*w_e)**2)*mu(x,y)
D = (4*(B_e**3))/(w_e**2)
#Calculating vibrational frequency for other isotopes.
def isotope_effect_rule(w_e_1, mu_1, mu_2):
    '''returns w_e_2*'''
    return (w_e_1)/(mu_1/mu_2)**.5

'''Statistical Mechanics Portion'''
#constants
B = 10.59341
h = 6.626*10**(-34) #J*sec
c = c_speed
k_B = 1.3807*10**(-23) #J*sec / K

def N_J(J, T):
    '''returns the number of molecules in a given rotational state
    at a given temperature'''
    return (2*J +1)*(e**(-(B*J*(J+1)*h*c)/(k_B*T)))
#RT
# print("RT_1 is", N_J(1, 300))
# print("100_2 is", N_J(2, 300))
# print("RT_3 is", N_J(3, 300))
# print("RT_4 is", N_J(4, 300))
# print("RT_5 is", N_J(5, 300))
# print("RT_6 is", N_J(6, 300))
# print("RT_7 is", N_J(7, 300))
# print("RT_8 is", N_J(8, 300))
# print("RT_9 is", N_J(9, 300))
# print("RT_10 is", N_J(10, 300))

#100 K
# print(N_J(1, 100))
# print(N_J(2, 100))
# print(N_J(3, 100))
# print(N_J(4, 100))
# print(N_J(5, 100))
# print(N_J(6, 100))
# print(N_J(7, 100))
# print(N_J(8, 100))
# print(N_J(9, 100))
# print(N_J(10, 100))

#500K
# print(N_J(1, 500))
# print(N_J(2, 500))
# print(N_J(3, 500))
# print(N_J(4, 500))
# print(N_J(5, 500))
# print(N_J(6, 500))
# print(N_J(7, 500))
# print(N_J(8, 500))
# print(N_J(9, 500))
# print(N_J(10, 500))

#1000 K
print(N_J(1, 1000))
print(N_J(2, 1000))
print(N_J(3, 1000))
print(N_J(4, 1000))
print(N_J(5, 1000))
print(N_J(6, 1000))
print(N_J(7, 1000))
print(N_J(8, 1000))
print(N_J(9, 1000))
print(N_J(10, 1000))
