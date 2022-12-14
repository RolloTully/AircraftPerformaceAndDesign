
from numpy import e, array, float32, int_, longdouble
import numpy as np
import matplotlib.pyplot as plt
import csv
from tkinter import filedialog

class Atmosphere(object):
    def __init__(self, g = 9.8065, R = 287.05287, Gamma = 1.4, L0 = -0.0065 ,L11 = 0, L20 = 0.001, T0 =288.15, T11 = 216.65, T20 = 216.65, P0 =101325, P11 = 22632.559, RF = 0.1): #Defaults to ISA
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
        self.Rho0= 1.225
        self.RF = RF #fuel reserve fraction

    def get_AtmosProperties(self, height):
        '''Computes Static pressure, Temperature, density, and the local speed of sound'''
        if 0.0 <= height <= 11e3:
            self.Local_Temperature = self.T0+self.L0*(height-0)
        elif 11e3 < height <= 20e3:
            self.Local_Temperature = self.T11+self.L11*(height-11)
        elif 20e3 < height <= 32e3:
            self.Local_Temperature = self.T20+self.L20*(height-20)
        elif 32e3< height:
            print("Height out of range")
        if 0.0 <= height <= 11e3:
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
                    self.formatted = np.array(self.formatted).astype(longdouble)
                    self.ExportFlightPlan.append(self.formatted)
            return array(self.ExportFlightPlan)
    def Mach_from_pressures(self, Sp, Ip,gamma):
        return  ((2/(gamma-1))*(((Ip/Sp)+1)**((gamma-1)/gamma)-1))**0.5



class Craft(object):
    def __init__(self, Reference_Area, k, cd0, mtow, fmf, n2, ct2):
        self.atmosphere = Atmosphere()
        self.tools = Tools()




        self.FMF = fmf

        self.TSFCLE = n2
        self.TSFCLC = ct2
        self.Fuel_Density = 785
        self.passenger_number = 2
        self.OEW = 134000*9.81
        self.MTOW = 227000*9.81
        self.Max_Payload = 44000*9.81
        self.Max_Fuel = 73300*9.81
        self.Cruise = 0.7
        self.Static_Thrust = 2*2570000
        self.S = 260
        self.Cd0 = 0.3
        self.K = 0.291
        self.SFC = 0.0165
        self.max_cruise_alt = 31000

        self.FW = 0.8*self.Max_Fuel
        self.DW = self.OEW+self.Max_Payload
        self.V_Ne = 203
        self.M_Ne = 0.88
        self.Tk_Cl = 2.44
        self.L_Cl = 2.98
        self.q_limit = 101325*(( 1+ (0.4/2)*(self.V_Ne/331)**2   )**(1.4/(0.4)) -1)
    def Block_Fuel_Range(self, m, h):
        self.cl_max_sar = (self.Cd0/(3*self.K))**0.5
        self.block_fuel_points = np.linspace(0,self.FMF, 100)
        self.Alt_Block_Fuel_Range = [self.Breguet_Altitude(self.cl_max_sar, m, h, self.ratio) for self.ratio in self.block_fuel_points]
        self.Cl_Block_Fuel_Range = [self.Breguet_Mach(self.cl_max_sar, m, h, self.ratio) for self.ratio in self.block_fuel_points]
        self.Mach_Block_Fuel_Range = [self.Breguet_Cl(self.cl_max_sar, m, h, self.ratio) for self.ratio in self.block_fuel_points]
        return self.Alt_Block_Fuel_Range, self.Cl_Block_Fuel_Range, self.Mach_Block_Fuel_Range


    def Breguet_Altitude(self, cl, m,h,zeta, operating_weight): # Working
        #print(zeta,operating_weight)
        self.w_e = operating_weight * (1-zeta)
        self.w_i = operating_weight
        self.atmo = self.atmosphere.get_AtmosProperties(h)
        self.clcd = (cl**0.5)/(self.Cd0 + self.K*cl**2)
        #print(self.w_i, self.w_e, self.w_i/self.w_e, zeta)
        return ((m*self.atmo[3])/9.81)*(1/self.SFC)*(cl/(self.Cd0+self.K*cl**2))*np.log(self.w_i/self.w_e)


    def Breguet_Mach(self, cl, m, h, zeta, operating_weight): #Working
        self.w_e = operating_weight * (1-zeta)
        self.w_i = operating_weight
        self.atmo = self.atmosphere.get_AtmosProperties(h)
        self.clcd = cl/(self.Cd0 + self.K*cl**2)
        return (2/(9.81*self.SFC))*((2*self.w_i)/(self.atmo[2]*self.S*cl))**0.5*self.clcd*(1-(self.w_i/self.w_e)**-0.5)

    def Breguet_Cl(self, cl, m, h, zeta, operating_weight): #Working
        self.W_e = operating_weight * (1-zeta)
        self.W_i = operating_weight
        self.atmo = self.atmosphere.get_AtmosProperties(h)
        self.clcd = cl**0.5/(self.Cd0 + self.K*cl**2)
        return ((self.atmo[3]*m)/(9.81*self.SFC))*(1/(self.K*self.Cd0))**0.5*(np.arctan((self.W_i/(0.5*self.atmo[2]*(m*self.atmo[3])**2*self.S))*(self.K/self.Cd0)**0.5)-np.arctan(self.W_e/(0.5*self.atmo[2]*(m*self.atmo[3])**2*self.S))*(self.K/self.Cd0)**0.5)

    def SAR_to_mpg(self, SAR):
        return ((((SAR/1609)*840))/4.546)*self.passenger_number #miles per kg

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
        elif gamma<0:
            self.thrust_req = self.drag-(self.FW+self.DW)*9.81*np.sin(gamma)
        return self.thrust_req

    def TSFC(self, t, m):#working
        return (1/(3600))*(self.TSFCLC*np.sqrt(t/self.atmosphere.T0)*m**self.TSFCLE)

    def Flight_Envelope(self):
        self.alt_array = range(1,40000)
        self.Mach_limit = []
        self.q_max = []
        self.Stall_limit_Tk = []
        self.Stall_limit_L = []
        self.V_max_H = []
        self.V_max_T = []
        self.t_a=[]
        #self.V_max_mach_limit = [self.cruise*self.atmosphere.get_AtmosProperties(h*0.3048) for h in self.alt_array]
        for self.ft_alt in range(1,40000):
            self.alt = self.ft_alt*0.3048
            self.atmo = self.atmosphere.get_AtmosProperties(self.alt)
            self.thrust_avalable = self.Static_Thrust*(self.atmo[2]/self.atmosphere.Rho0)
            self.Mach_limit.append(self.M_Ne*np.sqrt(1.4*287*self.atmo[0]))
            self.Stall_limit_Tk.append(np.sqrt((2/self.atmo[2])*(self.MTOW/self.S)*(1/self.Tk_Cl)))
            self.Stall_limit_L.append(np.sqrt((2/self.atmo[2])*(self.MTOW/self.S)*(1/self.L_Cl)))
            self.q_max.append(np.sqrt((2*self.q_limit)/self.atmo[2]))
            print(np.sqrt((2*self.q_limit)/self.atmo[2]))
            print(self.q_limit, self.atmo[2])
            #self.Stall_limit_L.append(np.sqrt((2/self.atmo[2])*(self.MTOW/self.S)*(1/self.L_Cl)))
            #self.V_max_H.append( ((1/(self.atmo[3]*self.Cd0))*((self.Static_Thrust/(self.MTOW))*( (self.MTOW)/self.S)+( (self.MTOW)/self.S)*np.sqrt( (((self.Static_Thrust/(self.MTOW)))**2)-4*self.Cd0*self.K  )   )  )**0.5)
            #self.V_max_T.append( ((1/(self.atmo[3]*self.Cd0))*((self.Static_Thrust/(self.MTOW))*( (self.MTOW)/self.S)-( (self.MTOW)/self.S)*np.sqrt( (((self.Static_Thrust/(self.MTOW)))**2)-4*self.Cd0*self.K  )   )  )**0.5)
            self.V_max_T.append( ((1/(self.atmo[2]*self.Cd0))*((self.thrust_avalable/(self.MTOW))*( (self.MTOW)/self.S)+( (self.MTOW)/self.S)*np.sqrt((((self.thrust_avalable/(self.MTOW)))**2)-4*self.Cd0*self.K  )   )  )**0.5)
            self.V_max_H.append( ((1/(self.atmo[2]*self.Cd0))*((self.thrust_avalable/(self.MTOW))*( (self.MTOW)/self.S)-( (self.MTOW)/self.S)*np.sqrt((((self.thrust_avalable/(self.MTOW)))**2)-4*self.Cd0*self.K  )   )  )**0.5)
            # ((1/(self.atmo[3]*self.Cd0))*((self.Static_Thrust/(self.MTOW))*( (self.MTOW)/self.S)+( (self.MTOW)/self.S)*np.sqrt( (((self.Static_Thrust/(self.MTOW)))**2)-4*self.Cd0*self.K  )   )  )**0.5
        '''Formatting Plots'''
        #plt.grid()
        #plt.plot(self.t_a)
        #plt.show()
        plt.title("Calculated flight evnelope for A300-700")
        plt.plot(self.Mach_limit,self.alt_array,label = "Mach Limit")
        plt.plot(self.Stall_limit_Tk ,self.alt_array,linestyle = '--',  label = "Takeoff Stall Limit")
        plt.plot(self.Stall_limit_L ,self.alt_array,linestyle = '--',  label = "Landing Stall Limit")
        plt.plot(self.V_max_T, self.alt_array, color = "r", label = "Max thrust available speed")
        plt.plot(self.V_max_H, self.alt_array, color = "r")
        plt.plot(self.q_max,self.alt_array, color = 'violet', label = "Structural limit")
        plt.ylim([0,None])
        plt.legend()
        plt.text(30,self.max_cruise_alt+200,"Max cruise altitude")
        plt.axhline(self.max_cruise_alt,color = 'k', linestyle = '-.', label = "Max Cruise altitude")
        plt.show()





class main():
    def __init__(self):
        self.atmosphere = Atmosphere()
        self.craft = Craft(361.6, 0.0259, 0.0417, 78000, 0.545, 0.432, 0.611)
        self.tools = Tools()

        self.FlightPlan = self.tools.ImportFlightPlan()
        self.Time = self.FlightPlan[:,0]#-self.FlightPlan[0,0]
        self.Time = self.Time# - self.Time[0]
        self.Latitude = self.FlightPlan[:,1]
        self.Longditude = self.FlightPlan[:,2]
        self.Alt  = self.FlightPlan[:,3]/3.280839895 #Feet
        #self.ImpactPressure = self.FlightPlan[:,4]
        self.speed = self.FlightPlan[:,4]*0.514
        #self.Temperature = self.FlightPlan[:,5]
        self.craft.Flight_Envelope()
        self.mainloop()
        self.Payload_Range_Chart()

    def Payload_Range_Chart(self):
        self.cl = (self.craft.Cd0/(3*self.craft.K))**0.5
        self.Ferry_Zeta = self.craft.Max_Fuel/(self.craft.OEW+self.craft.Max_Fuel)

        self.Payload_Zeta = self.craft.Max_Fuel/(self.craft.OEW+self.craft.Max_Payload+self.craft.Max_Fuel)
        self.OEWPAY_weight = []
        self.OEWPAYRES_weight = []
        self.Total_weight = []
        self.range = []
        '''First section, increasing fuel'''
        self.block_fuel_ratio_range  = [0,(self.craft.OEW+self.craft.Max_Payload)/self.craft.MTOW]
        for self.Fuel_weight in np.linspace(0,self.craft.MTOW-(self.craft.OEW+self.craft.Max_Payload),100):
            self.OEWPAY_weight.append(self.craft.OEW+self.craft.Max_Payload)
            self.OEWPAYRES_weight.append(self.craft.OEW+self.craft.Max_Payload+0.1*self.Fuel_weight)
            self.Total_weight.append(self.craft.OEW+self.craft.Max_Payload+self.Fuel_weight)
            #print(self.craft.MTOW,self.craft.OEW,self.craft.Max_Payload)
            self.range.append(self.craft.Breguet_Altitude(self.cl,self.craft.Cruise,10668,(self.Fuel_weight/(self.craft.OEW+self.craft.Max_Payload+self.Fuel_weight)),self.Fuel_weight+self.craft.OEW+self.craft.Max_Payload))
        plt.axvline(self.craft.Breguet_Altitude(self.cl,self.craft.Cruise,10668,(self.Fuel_weight/(self.craft.OEW+self.craft.Max_Payload+self.Fuel_weight)),self.Fuel_weight+self.craft.OEW+self.craft.Max_Payload), color = 'g', linestyle = '--')
        plt.text(self.craft.Breguet_Altitude(self.cl,self.craft.Cruise,10668,(self.Fuel_weight/(self.craft.OEW+self.craft.Max_Payload+self.Fuel_weight)),self.Fuel_weight+self.craft.OEW+self.craft.Max_Payload), 5000, "Max Payload Range", rotation = 270)
        print("Max Economic Range:")
        print(self.craft.Breguet_Altitude(self.cl,self.craft.Cruise,10668,(self.Fuel_weight/(self.craft.OEW+self.craft.Max_Payload+self.Fuel_weight)),self.Fuel_weight+self.craft.OEW+self.craft.Max_Payload))
        '''Second section, decreasing payload, increasing fuel'''
        self.F_Fuel_weight = self.Fuel_weight
        for self.Fuel_weight in np.linspace(self.F_Fuel_weight, self.craft.Max_Fuel, 100):
            self.OEWPAYRES_weight.append(self.craft.OEW+(self.craft.MTOW-(self.craft.OEW+self.Fuel_weight*0.9)))
            self.OEWPAY_weight.append(self.craft.OEW+(self.craft.MTOW-(self.craft.OEW+self.Fuel_weight)))
            self.Total_weight.append(self.craft.MTOW)
            self.zeta = self.Fuel_weight/self.craft.MTOW
            self.range.append(self.craft.Breguet_Altitude(self.cl,self.craft.Cruise,10668,self.zeta,self.craft.MTOW))
        plt.axvline(self.craft.Breguet_Altitude(self.cl,self.craft.Cruise,10668,self.zeta,self.craft.MTOW), color = 'b', linestyle = '--')
        plt.text(self.craft.Breguet_Altitude(self.cl,self.craft.Cruise,10668,self.zeta,self.craft.MTOW), 5000, "Max Economic Range", rotation = 270)
        print("Max Payload range:")
        print(self.craft.Breguet_Altitude(self.cl,self.craft.Cruise,10668,self.zeta,self.craft.MTOW))

        '''Third section decreasing payload, constant fuel'''
        self.starting_payload_weight = self.craft.MTOW-(self.craft.OEW+self.craft.Max_Fuel)
        for self.Payload_weight in np.linspace(self.starting_payload_weight,0,100):
            self.OEWPAYRES_weight.append(self.craft.OEW+self.Payload_weight+0.1*self.craft.Max_Fuel)
            self.OEWPAY_weight.append(self.craft.OEW+self.Payload_weight)
            self.Total_weight.append(self.craft.OEW+self.craft.Max_Fuel+self.Payload_weight)
            self.zeta = self.craft.Max_Fuel/(self.craft.OEW+self.craft.Max_Fuel+self.Payload_weight)
            self.range.append(self.craft.Breguet_Altitude(self.cl,self.craft.Cruise,10668,self.zeta,self.craft.OEW+self.craft.Max_Fuel+self.Payload_weight))
        print("Maximum ferry range:")
        print(self.craft.Breguet_Altitude(self.cl,self.craft.Cruise,10668,self.zeta,self.craft.OEW+self.craft.Max_Fuel+self.Payload_weight))


        '''Plotting'''
        plt.title("Payload-Range chart for B777-200")
        plt.xlabel("Trip Range, (km)")
        plt.ylabel("Weight W,(kgf)")
        plt.grid()

        plt.plot(self.range,self.OEWPAY_weight, label = "OEW + PAY")
        plt.plot(self.range,self.OEWPAYRES_weight, label = "OEW +PAY + RES")
        plt.plot(self.range,self.Total_weight, label = "OEW + PAY + FUEL")
        plt.fill_between(np.array(self.range)[::-1],np.array(self.Total_weight)[::-1],np.array(self.OEWPAYRES_weight)[::-1],np.array(self.Total_weight)[::-1]>np.array(self.OEWPAYRES_weight)[::-1] ,color = 'green', alpha = 0.15)
        plt.fill_between(np.array(self.range)[::-1],np.array(self.OEWPAYRES_weight)[::-1],np.array(self.OEWPAY_weight)[::-1],np.array(self.OEWPAYRES_weight)[::-1]>np.array(self.OEWPAY_weight)[::-1] ,color = 'orange', alpha = 0.15)
        plt.fill_between(np.array(self.range)[::-1],np.array(self.OEWPAY_weight)[::-1],np.full((300),self.craft.OEW)[::-1],np.array(self.OEWPAY_weight)[::-1]>np.full((300),self.craft.OEW)[::-1] ,color = 'blue', alpha = 0.15)
        plt.fill_between(np.array(self.range)[::-1],0,np.full((300),self.craft.OEW)[::-1],color = 'black', alpha = 0.15)
        plt.axvline(x=self.craft.Breguet_Altitude(self.cl, self.craft.Cruise, 10668,self.Ferry_Zeta, self.craft.OEW+self.craft.Max_Fuel), color = 'r', linestyle = '--')
        plt.text(self.craft.Breguet_Altitude(self.cl, self.craft.Cruise, 10668,self.Ferry_Zeta, self.craft.OEW+self.craft.Max_Fuel), 5000, "Max Ferry Range", rotation = 270)

        plt.axhline(y = self.craft.OEW, color = 'k', linestyle = '-.',label = "Operating empty weight")
        plt.axhline(y = self.craft.MTOW, color = 'k', linestyle = '-.', label = "Maximum takeoff weight")
        plt.text(100, self.craft.OEW+1500, "OEW")
        plt.text(100, self.craft.MTOW+1500, " MTOW")

        plt.ylim([0,self.craft.MTOW*1.1])
        plt.xlim([0,None])
        #plt.legend()
        plt.show()



        #print(self.Ferry_Zeta, self.Economic_Zeta, self.Payload_Zeta)
        #print()
        #print(self.craft.Breguet_Altitude(self.cl, 0.8, 10668,self.Economic_Zeta))
        #print(self.craft.Breguet_Altitude(self.cl, 0.8, 10668,self.Payload_Zeta))





    def mainloop(self):
        self.Atmos_conditions = np.array([self.atmosphere.get_AtmosProperties(alt) for alt in self.Alt]) #Computes the atmospheric properties alon the flight path
        #self.mach = np.array([self.tools.Mach_from_pressures(self.slice[0], self.slice[1], self.atmosphere.gamma) for self.slice in np.dstack((self.Atmos_conditions[:,1], self.ImpactPressure))[0]])
        self.mach = self.speed/self.Atmos_conditions[:,3]
        self.TAS = self.speed # We make this assumtion as there is no pitotstatic data avalable
        #self.TAS = self.mach*self.Atmos_conditions[:,3]
        self.climb_rate = (np.diff(self.Alt)/np.diff(self.Time)) #Computes the aircraft climb rate
        #print(self.climb_rate)
        #print(self.TAS)

        self.gamma = np.arcsin(self.climb_rate/self.TAS[1:])  #Flight path angle

        self.clcd_contin = np.array([self.craft.get_trimCLCD(self.slice[0],self.slice[1],self.slice[2]) for self.slice in np.dstack((self.gamma, self.mach[1:], self.Alt[1:]))[0]])

        self.t_req = np.array([self.craft.thurst_required(self.slice[0],self.slice[1],self.slice[2],self.slice[3]) for self.slice in np.dstack((self.gamma, self.Alt[1:], self.clcd_contin[:,1],self.mach[1:]))[0]])

        self.TSFC = np.array([self.craft.TSFC(self.slice[0], self.slice[1])  for self.slice in np.dstack((self.Atmos_conditions[:,0],self.mach))[0]])
        self.fuel_flow = (self.t_req/9.81)*self.TSFC[1:]
        self.fuel_usage = self.fuel_flow*np.diff(self.Time)
        self.comp_fw = np.array([self.craft.DW+self.craft.FW- np.sum(self.fuel_usage[0:index])   for index in range(0,self.fuel_usage.shape[0])])
        self.SAR = np.array([self.craft.get_SAR(self.slice[0], self.slice[1], self.slice[2]) for self.slice in np.dstack((self.gamma,self.mach[1:], self.Alt[1:]))[0]])
        self.SE = np.array([self.craft.get_SE(self.slice[0], self.slice[1], self.slice[2]) for self.slice in np.dstack((self.gamma,self.mach[1:], self.Alt[1:]))[0]])
        self.mpg = np.array([self.craft.SAR_to_mpg(sar) for sar in self.SAR])
        self.cl_max_sar = (self.craft.Cd0/(3*self.craft.K))**0.5
        #plt.plot(self.Time,self.mach)
        #plt.plot(self.Time,self.Alt/np.max(self.Alt))
        #plt.show()
        plt.plot(self.Time[1:], self.SAR)
        #plt.plot(self.SE)
        #plt.show()
        #plt.plot(self.Alt)
        #plt.show()
        #plt.plot(self.climb_rate)
        plt.title("Calculated SAR for B777-200")
        plt.xlabel("Travel time, Seconds")
        plt.ylabel("SAR, kg per meter")
        plt.grid()
        plt.show()

        '''
        self.Br_alt, self.Br_cl, self.Br_Mach = self.craft.Block_Fuel_Range(self.craft.Cruise, 10972)
        plt.plot(self.Br_alt)
        #plt.show()
        plt.plot(self.Br_cl)
        #plt.show()
        plt.plot(self.Br_Mach)
        plt.show()
        '''
        #self.Payload_Range_Chart()

if __name__=="__main__":
    main()
