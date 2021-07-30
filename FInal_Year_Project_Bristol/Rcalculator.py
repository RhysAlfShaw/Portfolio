# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 20:15:17 2019

@author: Rhys Alfred Shaw

This code was created for use in my Final Year Project, Where it would calculcuate the size of the galaxy clusters radius
R500 given kT(clusters ICM temperature) and its redshift. 
Its a GUI made with tkinter to help the work flow when analyseing cluster observations. 

On the most basic level its a custom GUI calculator.

"""
#import cosmolopy as cos
from astropy.cosmology import WMAP9 as cosmo
import numpy as np
import tkinter as tk
from tkinter import ttk , Entry, Label, Tk, Button

G = 4.30091E-3 #pc*Mo**(-1)(km/s)**2
#k = 8.617333262145E−5 #ev*K**-1
#kT = 14.4626 #K change for each iteration
#z = 0.596 #needs to change for each cluster
h = 0.696   
ho = 69.6   
alpha = 1.53 
mo = 3.02E14  
wm = 0.286  #omega matter
wv = 0.714  #omega delta        


class GUI:
    def __init__(self,master):
        self.master = master
        master.title("Calculator For R_500!")
        
        vcmd = master.register(self.validate)
        
        self.label = Label(master,text='welcome to this caluculator for R_500, Please enter values in the correct entries.')
        self.label.pack()
        
        self.zlabel = Label(master,text = 'Please enter the clusters redshift here:')
        self.zlabel.pack()
        self.z = Entry(master, validate='key', validatecommand=(vcmd,'%P'))
        self.z.pack()
        
        self.kTlabel = Label(master, text='Please enter the Clusters Temperture in (KeV):')
        self.kTlabel.pack()
        
        self.kT = Entry(master, validate='key',validatecommand=(vcmd,'%P'))
        self.kT.pack()
        
        self.calc = Button(master, text='---calculate---', command=lambda:self.command(master,'calculate'))
        self.calc.pack()
        
        self.clear = Button(master, text='--clear--', command=lambda:self.command(master,'clear'))
        self.clear.pack()
     
    def command(self,master,method):
        self.master = master
        if method =='calculate':
            
            kT = float(self.kT.get())   #takes kT (temperature in units of energy) of cluster
            z = float(self.z.get())     #takes redshift of cluster
            scale = cosmo.kpc_proper_per_arcmin(z)/60*10**(-3)
            
            hz = ho*(0.3*(1+z)**3 + 0.7)**0.5
            #hz = ho*(1+z)**(3/2)
            E = hz/ho
            
            M_500 = mo*(kT/5)**alpha * E**(-1)*(h)**(-1)

            pc_z = 3*(hz)**2/(8*np.pi*G*1E-6)
            
            R_500 = (3*M_500/(4*np.pi*500*pc_z))**(1/3)  #*(h)**(-1/3)
            
            arcsec = R_500/scale

            pixelrad = arcsec/0.492
            
            ld = cosmo.luminosity_distance(z)
            
            self.CD = Label(master,text=('Critical density at z:   {0:.3e}'.format(pc_z),'Mo/Mpc^3'))
            self.CD.pack()
            self.CM = Label(master,text=('Cluster Mass within R_500:   {0:.3e}'.format(M_500),'Mo'))
            self.CM.pack()
            self.R = Label(master,text=('R_500: {0:.3e}'.format(R_500),'Mpc'))
            self.R.pack()
            self.a = Label(master,text=(arcsec,'Arcsec'))
            self.a.pack()
            self.PR = Label(master, text=('Pixel radius : {0:.3e}'.format(pixelrad)))
            self.PR.pack()
            
            self.ld = Label(master, text=('Luminosity Distance : {0:.3e}'.format(ld)))
            self.ld.pack()
        if method =='clear':
            self.CD.pack_forget()
            self.CM.pack_forget()
            self.R.pack_forget()
            self.a.pack_forget()
            self.PR.pack_forget()
            self.ld.pack_forget()
            
            
           
    
    def validate(self, new_text):
        if not new_text: # the field is being cleared
             self.entered_number = 0
             return True
        try:
             self.entered_number = float(new_text)
             return True
        except ValueError:
             return False

root = Tk()
one = GUI(root)
root.mainloop()
