
from numpy import e, array, float32
import numpy as np
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
        if 0 < height <= 11e3:
            self.Local_Temperature = self.T0+self.L0*(height-0)
        elif 11e3 < height <= 20e3:
            self.Local_Temperature = self.T11+self.L11*(height-11)
        elif 20e3 < height <= 32e3:
            self.Local_Temperature = self.T20+self.L20*(height-20)
        elif 32e3< height:
            print("Height out of range")
        if 0 < height <= 11e3:
            self.Local_Pressure = self.P0*(1+(self.L0/self.T0)*height)**(-self.g/(self.R*self.L0))
        elif 11e3 < height <= 20e3:
            self.Local_Pressure = self.P11*e**((-self.g/(self.R*self.T11))*(height-11e3))
        elif 20e3 < height <= 32e3:
            '''To be Implemented'''
            pass
        elif 32e3 < height:
            print("Height out of range")
        self.Local_Desnity = self.Local_Pressure/(self.R*self.Local_Temperature)
        self.a = (self.gamma*self.R*self.Local_Temperature)**0.5
        return self.Local_Temperature, self.Local_Pressure, self.Local_Desnity, self.a

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

    def Mach_from_pressures(self, param,gamma):
        Sp,Ip = param[0],param[1]
        return  ((2/(gamma-1))*(((Ip/Sp)+1)**((gamma-1)/gamma)-1))**0.5

    def ImpactPressure(self, param, Gamma):
        return param[0]*((1+((Gamma-1)/2)*param[1]**2   )**(Gamma/(Gamma-1))-1)



class main():
    def __init__(self):
        self.atmosphere = Atmosphere()
        self.tools = Tools()
        self.FlightPlan = self.tools.ImportFlightPlan()
        self.Time = self.FlightPlan[:,0]
        self.Alt  = self.FlightPlan[:,3]
        self.ImpactPressure = self.FlightPlan[:,4]
        self.Temperature = self.FlightPlan[:,5]
        self.mainloop()
    def mainloop(self):
        self.Atmospressure = np.array([self.atmosphere.get_AtmosProperties(alt) for alt in self.Alt])
        #Static Pressure, Impact pressure
        self.param_pairs = np.dstack((self.Atmospressure[:,1], self.ImpactPressure))
        self.mach = np.array([self.tools.Mach_from_pressures(self.pair, self.atmosphere.gamma) for self.pair in self.param_pairs[0]])

        plt.plot(self.Time, self.Atmospressure[:,0]/np.amax(self.Atmospressure[:,0]))
        plt.plot(self.Time, self.Atmospressure[:,1]/np.amax(self.Atmospressure[:,1]))
        plt.plot(self.Time, self.Atmospressure[:,2]/np.amax(self.Atmospressure[:,2]))
        plt.plot(self.Time, self.Atmospressure[:,3]/np.amax(self.Atmospressure[:,3]))
        plt.plot(self.Time, self.Temperature/np.amax(self.Temperature))
        plt.plot(self.Time, self.mach)
        plt.show()

        plt.plot(self.Time, self.mach*self.Atmospressure[:,3])
        plt.plot(self.Time, self.mach*(self.Temperature*self.atmosphere.gamma*self.atmosphere.R)**0.5)
        plt.show()
if __name__=="__main__":
    main()
