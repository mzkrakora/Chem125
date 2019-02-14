''' Name = Marty Krakora
Date = 13 February 2019
Title = HCl DCl analysis '''

from math import *

# After fitting data points to curve, fill in values here to get the constants of B_0, B_1, v_1, etc.
#Note here we are detecting the overtone transition so m = 2
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


print("For H_35_Cl, v_0 = {v_0}, v_1 = {v_1}, alpha_e = {alpha_e}, B_0 = {B_0}, B_1 = {B_1}, r_0 = {r_0*100} ang, r_1 = {r_1*100} Ang, r_e = {r_e*10**10}, B_e = {B_e}, x_e = {x_e}, w_e = {w_e} (cm-1), k = {k}, D = {D}")
