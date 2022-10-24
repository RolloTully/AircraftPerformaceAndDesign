
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
            assert "Height out of range"
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
    def Mach_from_pressures(self, Sp, Ip,gamma):
        return  ((2/(gamma-1))*(((Ip/Sp)+1)**((gamma-1)/gamma)-1))**0.5
class Craft(object):
    def __init__(self, Reference_Area, k, cd0, mtow, fmf, n2, ct2):
        self.atmosphere = Atmosphere()
        self.tools = Tools()
        self.S = Reference_Area
        self.K = k
        self.Cd0 = cd0
        self.MTOW = mtow
        self.FMF = fmf
        self.FW = self.MTOW*self.FMF
        self.DW = self.MTOW - self.FW
        self.TSFCLE = n2
        self.TSFCLC = ct2

    def get_live_weight(self,flightplan):
        self.history = []
        self.FlightPlan = flightplan
        for self.index in range(0, self.FlightPlan.shape[0]-1):
            self.local_atmos = self.atmosphere.get_AtmosProperties(self.FlightPlan[self.index,3])
            self.time_step = self.FlightPlan[self.index+1,0] - self.FlightPlan[self.index,0]
            self.mach = self.tools.Mach_from_pressures(self.local_atmos[1], self.FlightPlan[self.index, 4], 1.4)
            self.climb_rate = (self.FlightPlan[self.index+1,3] - self.FlightPlan[self.index,3])/self.time_step
            self.TAS = self.mach * self.local_atmos[3]
            self.gamma = np.arcsin(self.climb_rate/self.TAS)
            self.clcd = self.get_trimCLCD(self.gamma, self.mach, self.FlightPlan[self.index,3])
            self.t_req = self.thurst_required(self.gamma, self.FlightPlan[self.index,3],self.clcd[1],self.mach)
            self.tsfc = self.TSFC(self.local_atmos[0], self.mach)
            self.fuel_rate = (self.t_req)*self.tsfc*(1/3600)
            print(self.t_req, self.fuel_rate,self.tsfc)
            self.FW = self.FW-self.fuel_rate*self.time_step
            print(self.gamma, self.mach , self.FlightPlan[self.index,3])
            self.history.append(self.FW)
        return self.history

    def get_trimCLCD(self, gamma_fp, m, h):
        self.atmo = self.atmosphere.get_AtmosProperties(h) #Pressure,local Temprature, Density, Speed of sound
        self.q = 0.5*self.atmosphere.gamma*self.atmo[1]*m**2
        self.L = (self.FW+self.DW)*9.81*np.cos(gamma_fp)
        self.cl = self.L/(self.q*self.S)
        self.cd = self.Cd0 + self.K*self.cl**2
        return [self.cl,self.cd]

    def thurst_required(self, gamma, h, cd, m):
        self.atmo = self.atmosphere.get_AtmosProperties(h)  #Pressure,local Temprature, Density, Speed of sound
        self.drag = 0.5*self.atmo[2]*((self.atmo[3]*m)**2)*self.S*cd
        self.thrust_req = self.drag+(self.FW+self.DW)*9.81*np.sin(gamma)
        return self.thrust_req

    def TSFC(self, t, m):
        return (1/3600*9.81)*(self.TSFCLC*np.sqrt(t/self.atmosphere.T0)*m**self.TSFCLE)

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
        plt.plot(self.craft.get_live_weight(self.FlightPlan))
        plt.show()
        #self.mainloop()

    def mainloop(self):
        self.Atmos_conditions = np.array([self.atmosphere.get_AtmosProperties(alt) for alt in self.Alt]) #Computes the atmospheric properties alon the flight path
        self.mach = np.array([self.tools.Mach_from_pressures(self.slice[0], self.slice[1], self.atmosphere.gamma) for self.slice in np.dstack((self.Atmos_conditions[:,1], self.ImpactPressure))[0]])
        self.TAS = self.mach*self.Atmos_conditions[:,3]
        self.climb_rate = (np.diff(self.Alt)/np.diff(self.Time)) #Computes the aircraft climb rate
        self.gamma = np.arcsin(self.climb_rate/self.TAS[:-1])  #Flight path angle
        self.clcd_contin = np.array([self.craft.get_trimCLCD(self.slice[0],self.slice[1],self.slice[2]) for self.slice in np.dstack((self.gamma, self.mach[:-1], self.Alt[:-1]))[0]])
        self.t_req = np.array([self.craft.thurst_required(self.slice[0],self.slice[1],self.slice[2],self.slice[3]) for self.slice in np.dstack((self.gamma, self.Alt[:-1], self.clcd_contin[:,1],self.mach[:-1]))[0]])
        self.TSFC = np.array([self.craft.TSFC(self.slice[0], self.slice[1])  for self.slice in np.dstack((self.Temperature,self.mach))[0]])
        self.fuel_flow = (self.t_req/9.81)*self.TSFC[:-1]*4
        plt.plot(self.t_req)
        plt.plot(self.fuel_flow)
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
