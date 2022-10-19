
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
class Craft(object):
    def __init__(self, Reference_Area, k, cd0, mtow, fmf, n2, ct2):
        self.atmosphere = Atmosphere()
        self.S = Reference_Area
        self.K = k
        self.Cd0 = cd0
        self.MTOW = mtow
        self.FMF = fmf
        self.FW = self.MTOW*self.FMF
        self.DW = self.MTOW - self.FW
        self.TSFCLE = n2
        self.TSFCLC = ct2
    def get_trimCLCD(self, gamma_fp, m, h):
        self.atmo = self.atmosphere.get_AtmosProperties(h)
        self.cl = (2*(self.DW+self.FW)*np.cos(gamma_fp))/(self.atmo[2]*self.atmo[3]*m*self.S)
        self.cd = self.Cd0 + self.K*self.cl**2
        return [self.cl,self.cd]
    def thurst_required(self, gamma, h, cd, m):
        self.atmo = self.atmosphere.get_AtmosProperties(h)
        self.drag = 0.5*self.atmo[1]*self.atmo[3]*m*self.S*cd
        return (self.FW+self.DW)*np.cos(gamma)+self.drag
    def fuel_flow(self, t, m):
        return  (1/(3600*9.81))*self.TSFCLC*np.sqrt(t)*m**self.TSFCLE
class main():
    def __init__(self):
        self.atmosphere = Atmosphere()
        self.craft = Craft(363.10, 0.0259, 0.0221, 217000, 0.358, 0.432, 0.611)
        self.tools = Tools()
        self.FlightPlan = self.tools.ImportFlightPlan()
        self.Time = self.FlightPlan[:,0]
        self.Alt  = self.FlightPlan[:,3]
        self.ImpactPressure = self.FlightPlan[:,4]
        self.Temperature = self.FlightPlan[:,5]
        self.mainloop()

    def mainloop(self):
        self.Atmos_conditions = np.array([self.atmosphere.get_AtmosProperties(alt) for alt in self.Alt])
        self.Static_pressure = self.Atmos_conditions[:,1]
        self.a = self.Atmos_conditions[:,3]
        #Static Pressure, Impact pressure
        self.mach = np.array([self.tools.Mach_from_pressures(self.pair, self.atmosphere.gamma) for self.pair in np.dstack((self.Static_pressure, self.ImpactPressure))[0]])
        self.TAS = self.mach*self.a
        self.climb_rate = (np.diff(self.Alt)/np.diff(self.Time))*np.sign(np.diff(self.Alt))
        self.gamma = np.arcsin(self.climb_rate/self.TAS[:-1])  #Flight path angle
        self.clcd_contin = np.array([self.craft.get_trimCLCD(self.slice[0],self.slice[1],self.slice[2]) for self.slice in np.dstack((self.gamma, self.mach[:-1], self.Alt[:-1]))[0]])
        self.t_req = np.array([self.craft.thurst_required(self.slice[0],self.slice[1],self.slice[2],self.slice[3]) for self.slice in np.dstack((self.gamma, self.Alt[:-1], self.clcd_contin[:,1],self.mach[:-1]))[0]])
        print(self.t_req)
        self.TSFR = np.array([self.craft.fuel_flow(self.slice[0], self.slice[1])  for self.slice in np.dstack((self.Temperature,self.mach))[0]])
        self.fuel_flow = self.t_req*self.TSFR[:-1]

        plt.plot(self.t_req/np.max(self.t_req))
        plt.plot(self.fuel_flow/np.max(self.fuel_flow))
        plt.show()

        '''
        plt.plot(self.Time, self.Atmos_conditions[:,0]/np.amax(self.Atmos_conditions[:,0]))
        plt.plot(self.Time, self.Atmos_conditions[:,1]/np.amax(self.Atmos_conditions[:,1]))
        plt.plot(self.Time, self.Atmos_conditions[:,2]/np.amax(self.Atmos_conditions[:,2]))
        plt.plot(self.Time, self.Atmos_conditions[:,3]/np.amax(self.Atmos_conditions[:,3]))
        plt.plot(self.Time, self.Temperature/np.amax(self.Temperature))
        plt.plot(self.Time, self.mach)
        plt.show()
        plt.plot(self.Time, self.mach*self.Atmos_conditions[:,3])
        plt.plot(self.Time, self.mach*(self.Temperature*self.atmosphere.gamma*self.atmosphere.R)**0.5)
        plt.show()
        '''

if __name__=="__main__":
    main()
