
from numpy import e
import matplotlib.pyplot as plt
import csv
from tkinter import filedialog
class Atmosphere(object):
    def __init__(self, g = 9.8065, R = 287.05287, Gamma = 1.4, L0 = -0.0065 ,L11 = 0, L20 = 0.001, T0 =288.15, T11 = 216.65, T20 = 216.65, P0 =101325, P11 = 22632.559): #Defaults to ISA
        self.g = g    #ms^-2
        self.R = R #JKg^-1K^-1
        self.gamma = Gamma
        self.L0 = L0  #Km^-1
        self.L11 = L11     #Km^-1
        self.L20 =L20
        self.T0 = T0  #kelvins
        self.T11=T11     #kelvins
        self.P0=P0     #pascals
        self.P11=P11 #pascals
        self.Rho0=1.225
    def get_AtmosProperties(self, height):
        '''Computes Static pressure, Temperature, density, and the local speed of sound'''
        if 0 < height <= 11:
            self.Local_Temperature = self.T0+self.L0*(height-0)
        elif 11 < height <= 20:
            self.Local_Temperature = self.T11+self.L11*(height-11)
        elif 20 < height:
            self.Local_Temperature = self.T20+self.L20*(height-20)
        elif 32 < height:
            print("Height out of range")
        if 0 < height <= 11:
            self.Local_Pressure = self.P0*(1+(self.L0/self.T0)*height)**(-self.g/(self.R*self.L0))
        elif 11 < height <= 20:
            self.Local_Pressure = self.P11*e**((-self.g/(self.R*self.T11))*(height-11))
        elif 20 < height:
            '''To be Implemented'''
            pass
        elif 32 < height:
            print("Height out of range")
        self.Local_Desnity = self.Local_Pressure/(self.R*self.Local_Temperature)
        return self.Local_Temperature, self.Local_Pressure, self.Local_Desnity

class Tools():
    def ImportFlightPlan(self, filename):
            self.FlightPlanFile = filedialog.askopenfilename()
            with open(self.FlightPlanFile, newline='') as self.FlightPlanCSV:
                self.FlightPlan =csv.reader()

class main():
    def __init__(self):
        self.atmosphere = Atmosphere()
        self.mainloop()


    def mainloop(self):
        print(self.atmosphere.get_AtmosProperties(10))
if __name__=="__main__":
    main()
