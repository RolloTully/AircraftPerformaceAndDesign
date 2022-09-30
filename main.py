
from numpy import e, array, float32
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
        self.T20 = T20
        self.P0=P0     #pascals
        self.P11=P11 #pascals
        self.Rho0=1.225
    def get_AtmosProperties(self, height):
        '''Computes Static pressure, Temperature, density, and the local speed of sound'''
        if 0 < height <= 11000:
            self.Local_Temperature = self.T0+self.L0*(height-0)
        elif 11000 < height <= 20000:
            self.Local_Temperature = self.T11+self.L11*(height-11)
        elif 20000 < height <= 32000:
            self.Local_Temperature = self.T20+self.L20*(height-20)
        elif 32000 < height:
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
    def ImportFlightPlan(self):
            self.ExportFlightPlan = []
            self.FlightPlanFile = filedialog.askopenfilename()
            with open(self.FlightPlanFile, newline='') as self.FlightPlanCSV:
                self.FlightPlan = csv.reader(self.FlightPlanCSV, delimiter=' ', quotechar='|')
                for self.row in list(self.FlightPlan)[1:]:
                    self.formatted = self.row[0].split(",")
                    self.formatted = array(self.formatted, dtype = float32)
                    self.ExportFlightPlan.append(self.formatted)
            return array(self.ExportFlightPlan)
    def ImpactPressure(self, StaticPressure, Mach, Gamma):
        return StaticPressure*((1+((Gamma-1)/2)*Mach**2   )**(Gamma/(Gamma-1))-1)


class main():
    def __init__(self):
        self.atmosphere = Atmosphere()
        self.tools = Tools()
        self.FlightPlan = self.tools.ImportFlightPlan()
        print(self.FlightPlan)
        self.Time = self.FlightPlan[:,0]
        self.Alt  = self.FlightPlan[:,3]
        self.Mach = self.FlightPlan[0,4]
        self.Atmospressure = [self.atmosphere.get_AtmosProperties(alt) for alt in self.Alt]
        print(self.Atmospressure)
        self.ImpactPressures = self.tools.ImpactPressure(self.Atmospressure, self.Mach, self.gamma)
        #self.mainloop()
    def mainloop(self):
        print(self.atmosphere.get_AtmosProperties(10))
if __name__=="__main__":
    main()
