
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
        self.Fuel_Density = 0
        self.passenger_number = 100

    def SAR_to_mpg(self, SAR):
        return ((((SAR/1609)*840)/1000)/4.546)*self.passenger_number #miles per kg

    def get_SAR(self,gamma,m,h):
        self.atmo = self.atmosphere.get_AtmosProperties(h)
        self.v = m * self.atmo[3]
        self.cd_trim = self.get_trimCLCD(gamma,m,h)[1]
        return self.v/(self.thurst_required(gamma,h,self.cd_trim,m)*self.TSFC(self.atmo[0],m))

    def get_SE(self, gamma, m ,h):
        self.atmo = self.atmosphere.get_AtmosProperties(h)
        self.cd_trim = self.get_trimCLCD(gamma ,m ,h)[1]
        return 1/(self.thurst_required(gamma ,h ,self.cd_trim ,m)*self.TSFC(self.atmo[0] ,m))

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
        if gamma>=0:
            self.thrust_req = self.drag+(self.FW+self.DW)*9.81*np.sin(gamma)
        if gamma<0:
            self.thrust_req = self.drag-(self.FW+self.DW)*9.81*np.sin(gamma)
        return self.thrust_req

    def TSFC(self, t, m):#working
        return (1/(3600*9.81))*(self.TSFCLC*np.sqrt(t/self.atmosphere.T0)*m**self.TSFCLE)

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
        self.Atmos_conditions = np.array([self.atmosphere.get_AtmosProperties(alt) for alt in self.Alt]) #Computes the atmospheric properties alon the flight path
        self.mach = np.array([self.tools.Mach_from_pressures(self.slice[0], self.slice[1], self.atmosphere.gamma) for self.slice in np.dstack((self.Atmos_conditions[:,1], self.ImpactPressure))[0]])
        self.TAS = self.mach*self.Atmos_conditions[:,3]
        self.climb_rate = (np.diff(self.Alt)/np.diff(self.Time)) #Computes the aircraft climb rate
        self.gamma = np.arcsin(self.climb_rate/self.TAS[1:])  #Flight path angle
        self.clcd_contin = np.array([self.craft.get_trimCLCD(self.slice[0],self.slice[1],self.slice[2]) for self.slice in np.dstack((self.gamma, self.mach[1:], self.Alt[1:]))[0]])
        self.t_req = np.array([self.craft.thurst_required(self.slice[0],self.slice[1],self.slice[2],self.slice[3]) for self.slice in np.dstack((self.gamma, self.Alt[1:], self.clcd_contin[:,1],self.mach[1:]))[0]])
        self.TSFC = np.array([self.craft.TSFC(self.slice[0], self.slice[1])  for self.slice in np.dstack((self.Temperature,self.mach))[0]])
        self.fuel_flow = (self.t_req/9.81)*self.TSFC[1:]
        self.fuel_usage = self.fuel_flow*np.diff(self.Time)
        self.comp_fw = np.array([self.craft.DW+self.craft.FW- np.sum(self.fuel_usage[0:index])   for index in range(0,self.fuel_usage.shape[0])])
        self.SAR = np.array([self.craft.get_SAR(self.slice[0], self.slice[1], self.slice[2]) for self.slice in np.dstack((self.gamma,self.mach[1:], self.Alt[1:]))[0]])
        self.SE = np.array([self.craft.get_SE(self.slice[0], self.slice[1], self.slice[2]) for self.slice in np.dstack((self.gamma,self.mach[1:], self.Alt[1:]))[0]])
        self.mpg = np.array([self.craft.SAR_to_mpg(sar) for sar in self.SAR])
        #plt.plot(self.TSFC)
        #plt.show()
        #plt.plot(self.mpg)
        plt.plot(self.SAR)
        plt.show()
        plt.plot(self.SE)
        plt.show()
if __name__=="__main__":
    main()
