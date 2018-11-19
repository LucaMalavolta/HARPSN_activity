#! /usr/bin/python
# -*- coding: utf-8 -*-
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
#import seaborn as sns

__version__ = '2016-10-10'
__author__ = 'Luca Malavolta'

def line_flux(w_cent,w_size):
   cond = (np.abs(wave_rest-w_cent)<=w_size)
   wl_range = sum(dwave_rest[cond])
   cont = sum(data[cond]) / wl_range
   cont_sig  = sum((data[cond]))/ wl_range**2
   return cont, cont_sig

def line_integ(w_cent,w_size):
   cond = (np.abs(wave_rest-w_cent)<=w_size)
   wl_range = sum(dwave_rest[cond])
   cont = sum(data[cond]) / wl_range * w_size * 2.
   cont_sig  = sum((data[cond])) #* (w_size * 2. / wl_range)**2
   return cont, cont_sig

def highest10(w_cent,w_size):
   cond = (np.abs(wave_rest-w_cent)<=w_size)
   vals = data[cond]/dwave_rest[cond]
   ind = np.argpartition(vals, -10)[-10:]
   cont = np.average(vals[ind])
   cont_sig = (np.std(vals[ind])**2)/9.
   return cont, cont_sig

def line_flux_triangle(w_cent,w_size):
   cond = (np.abs(wave_rest-w_cent)<=w_size)
   triangle = 1. - np.abs(wave_rest-w_cent)/w_size
   denom = sum(triangle[cond]*dwave_rest[cond])
   print cond
   print w_cent
   print wave_rest
   print triangle[cond]
   print sum(triangle[cond]*dwave_rest[cond])
   print triangle
   cent = sum((data[cond])*triangle[cond])/denom
   cent_sig = sum(  (data[cond])*(triangle[cond])**2 ) /denom**2
   return cent, cent_sig

def add_subplot(plot_obj, title, wide='False', half_size=False):
   plot_obj.set_ylabel(title,fontsize=10)
   #plot_obj.title.set_text(title)

   # We change the fontsize of minor ticks label
   plot_obj.tick_params(axis='both', which='major', labelsize=10)
   plot_obj.tick_params(axis='both', which='minor', labelsize=8)

   if wide==True:
       print title, np.amin(wave_rest)
       xlim_min = np.amax([w_cent-w_size*5., np.amin(wave_rest)])
       xlim_max = np.amin([w_cent+w_size*5., np.amax(wave_rest)])
   else:
       xlim_min = w_cent-w_size*2.
       xlim_max = w_cent+w_size*2.
   sel = (wave_rest>xlim_min) & (wave_rest<xlim_max)
   #ylim = np.max(data[sel]/norm_factor/dwave_rest[sel])
   plot_obj.set_xlim(xlim_min,xlim_max)
   plot_obj.set_ylim(0.,1.5)
   plot_obj.axvline(w_cent-w_size/2.0, c='r')
   plot_obj.axvline(w_cent+w_size/2.0, c='r')
   if half_size:
      plot_obj.axvline(w_cent-w_size/4.0, c='g')
      plot_obj.axvline(w_cent+w_size/4.0, c='g')
   plot_obj.plot(wave_rest[sel],data[sel]/norm_factor/dwave_rest[sel])

if __name__ == "__main__":

    print 'DISCLAIMER: the value obtained for the NaID are not good for star-to-star '
    print 'comparison, since the continuum regions include several spectral lines '
    print 'A different approach than Gomes da Silva 2011+ was followed'
    print 'These value are meant only to assess the internal variation of the activity level'
    print 'of a star'


    parser = argparse.ArgumentParser(description='Generalized Lomb-Scargle periodogram.', add_help=False)
    argadd = parser.add_argument   # function short cut
    argadd('-?', '-h', '-help', '--help', help='show this help message and exit', action='help')
    argadd('input_file', nargs='?',help='Data file. If not specified example will be shown.')

    args = parser.parse_args().__dict__
    input_file = args.pop('input_file')

    file_rad, file_ext = os.path.splitext(input_file)

    plot_dir = './' + file_rad + '_plots'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)


    cc = 299792.458

    fileout = open(file_rad+'.output','w')
    fileout.write('file_name                    \tbjd           \trv      \trv_noise\tR_out   \tR_sig   \tHalpha  \tHa_sig  \tNaID    \tNaID_sig\tHeI     \tHeI_sig \tHalphaHF\tHaHF_sig\tCaIquiet\tCaIq_sig\n')
    fileout.write('#############################\t##############\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\n')
    #for ii_list in xrange(0,10):


    s1d_file = input_file


    s1d_hdu = fits.open(s1d_file)
    s1d_header = s1d_hdu[0].header
    data = s1d_hdu[0].data
    naxis1 = s1d_header['NAXIS1']

    cdelt1 = s1d_header['CDELT1']
    crval1 = s1d_header['CRVAL1']

    wave = np.arange(0,naxis1,1.00000,dtype=np.double)*cdelt1 + crval1

    wave_rest = wave.copy()
    print wave_rest
    dwave_rest = np.ones(naxis1) * cdelt1

    ## R H&K calculation
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)

    wave_ref = (np.abs(wave - 4000.00) < 50 )
    norm_factor = np.sum(data[wave_ref])/np.sum(dwave_rest[wave_ref])
    #plots are normalized acording to the flux of the chosen order

    w_cent = 3933.664
    w_size = 1.0900
    CaIIK, CaIIK_sig = line_flux_triangle(w_cent,w_size)
    add_subplot(ax1,'CaII K',True)

    w_cent = 3968.470
    w_size = 1.0900
    CaIIH, CaIIH_sig = line_flux_triangle(w_cent,w_size)
    add_subplot(ax2,'CaII H',True)

    w_cent = 3901.070
    w_size = 10.
    Vcont, Vcont_sig = line_flux(w_cent,w_size)
    add_subplot(ax3,'Left cont.')

    w_cent = 4001.070
    w_size = 10.
    Rcont, Rcont_sig = line_flux(w_cent,w_size)
    add_subplot(ax4,'Right cont.')

    R_raw=(CaIIK+CaIIH)/(Vcont+Rcont)
    R_out = 1.111 * R_raw + 0.0153
    R_sig = 1.111 * np.sqrt((CaIIK_sig+CaIIH_sig)/(Vcont+Rcont)**2 + (CaIIK+CaIIH)**2*(Vcont_sig+Rcont_sig)/(Vcont+Rcont)**4)

    print 'CaII K ', CaIIK, np.sqrt(CaIIK_sig)
    print 'CaII H ', CaIIH, np.sqrt(CaIIH_sig)
    print 'CaII V ', Vcont, np.sqrt(Vcont_sig)
    print 'CaII R ', Rcont, np.sqrt(Rcont_sig)
    print 'R out  ', R_out, R_sig
    print
    fig.savefig(plot_dir + '/CaII_' + input_file + '.png', dpi = 300)


    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig2, (ax4, ax5, ax6) = plt.subplots(3)

    wave_ref = (np.abs(wave - 6560.00) < 50 )
    norm_factor = np.sum(data[wave_ref])/np.sum(dwave_rest[wave_ref])
    #plots are normalized acording to the flux of the chosen order

    ## Halpha Central
    ## Central wavelength at 6562.828 (Gomez da Silva 2011+, Robertson 2013a+)
    w_cent = 6562.808
    w_size = 0.8
    order =  64
    Hcent, Hcent_sig = line_flux(w_cent,w_size)
    add_subplot(ax1,'Halpha cent',True,True)

    # first continuum window in the range [-700:-300] km/s, as in Robertson 2013a+, Kurster
    w_cent = 6550.887
    w_size = 5.375
    order =  64
    Lcont, Lcont_sig = line_flux(w_cent,w_size)
    add_subplot(ax2,'Left cont.',False)
    add_subplot(ax5,'Left cont.',False)

    w_cent = 6580.31
    w_size = 4.37825
    order =  64
    Rcont, Rcont_sig = line_flux(w_cent,w_size)
    add_subplot(ax3,'Right continuum',False)
    add_subplot(ax6,'Right continuum',False)

    Halpha =(Hcent)/(Lcont+Rcont)
    Halpha_sig = Halpha * np.sqrt(Hcent_sig/Hcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)

    print 'H cent ', Hcent, np.sqrt(Hcent_sig)
    print 'H Lref ', Lcont, np.sqrt(Lcont_sig)
    print 'H Rref ', Rcont, np.sqrt(Rcont_sig)
    print 'Halpha ', Halpha, Halpha_sig
    print
    fig.savefig(plot_dir + '/Halpha_' + input_file + '.png', dpi = 300)

    ## Halpha Central with half-size window
    ## Central wavelength at 6562.828 (Gomez da Silva 2011+, Robertson 2013a+)
    w_cent = 6562.808
    w_size = 0.4
    order =  64
    Hcent, Hcent_sig = line_flux(w_cent,w_size)

    HalphaHF =(Hcent)/(Lcont+Rcont)
    HalphaHF_sig = HalphaHF * np.sqrt(Hcent_sig/Hcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)

    print 'HalphaHF ', HalphaHF, HalphaHF_sig
    print



    ############### CaI activity-less index
    ## CaI "quiet" index (Kurster 2003+, Robertson 2013a+)


    w_cent = 6572.795
    w_size = 0.34
    order =  64
    CaIcent, CaIcent_sig = line_flux_triangle(w_cent,w_size)

    ## Rcont, Lcont as in the previous definition
    CaIquiet = 2.*(CaIcent)/(Lcont+Rcont)
    CaIquiet_sig = CaIquiet * np.sqrt(CaIcent_sig/CaIcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
    print 'CaI quiet ', CaIquiet, CaIquiet_sig
    print
    add_subplot(ax4,'CaI central',True)
    fig2.savefig(plot_dir + '/CaIquiet_' + input_file + '.png', dpi = 300)



    ############### NaI D alt
    ## Alternative approach to obtain Na I D1 & D2 lines
    ## check the disclaimer
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)
    wave_ref = (np.abs(wave - 5900.00) < 50 )
    norm_factor = np.sum(data[wave_ref])/np.sum(dwave_rest[wave_ref])
    #plots are normalized acording to the flux of the chosen order

    w_cent = 5895.920
    w_size = 0.25
    order =  53
    NaID1, NaID1_sig = line_flux(w_cent,w_size)
    #NaID1, NaID1_sig = line_integ(w_cent,w_size)
    add_subplot(ax1,'NaID1',True)

    w_cent = 5889.950
    w_size = 0.25
    order =  53
    NaID2, NaID2_sig = line_flux(w_cent,w_size)
    #NaID2, NaID2_sig = line_integ(w_cent,w_size)
    add_subplot(ax2,'NaID2',True)

    ## Reference region
    w_cent = 5805.000
    w_size = 5.0
    order =  52
    Lcont, Lcont_sig = line_flux(w_cent,w_size)
    #Lcont, Lcont_sig = line_integ(w_cent,w_size)
    add_subplot(ax3,'Left cont.')

    w_cent = 6090.000
    w_size = 10.00
    order =  57
    Rcont, Rcont_sig = line_flux(w_cent,w_size)
    #Rcont, Rcont_sig = line_integ(w_cent,w_size)
    add_subplot(ax4,'Right cont.')


    NaID =(NaID1+NaID2)/(Lcont+Rcont)
    NaID_sig = NaID * np.sqrt((NaID1_sig+NaID2_sig)/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)

    print 'NaID1 cen ', NaID1, np.sqrt(NaID1_sig)
    print 'NaID2 cen ', NaID2, np.sqrt(NaID2_sig)
    print 'NaID Lref ', Lcont, np.sqrt(Lcont_sig)
    print 'NaID Rref ', Rcont, np.sqrt(Rcont_sig)
    print 'NaID ', NaID, NaID_sig
    print
    fig.savefig(plot_dir + '/NaID_' + input_file + '.png', dpi = 300)

    ############### NaI D
    ## NaID computation
    ## as described in Sect. 3 of Gomes da Silva 2011+ and Dias 2007a+
    ## check the disclaimer

    w_cent = 5895.920
    w_size = 0.25
    order =  53
    NaID1, NaID1_sig = highest10(w_cent,w_size)

    w_cent = 5889.950
    w_size = 0.25
    order =  53
    NaID2, NaID2_sig = highest10(w_cent,w_size)

    ## Reference region

    w_cent = 5805.000
    w_size = 5.0
    order =  52
    Lcont, Lcont_sig = highest10(w_cent,w_size)

    w_cent = 6090.000
    w_size = 10.00
    order =  57
    Rcont, Rcont_sig = highest10(w_cent,w_size)


    NaID_alt =(NaID1+NaID2)/(Lcont+Rcont)
    NaID_alt_sig = NaID_alt * np.sqrt((NaID1_sig+NaID2_sig)/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)

    print 'NaID1_alt cen ', NaID1, np.sqrt(NaID1_sig)
    print 'NaID2_alt cen ', NaID2, np.sqrt(NaID2_sig)
    print 'NaID_alt Lref ', Lcont, np.sqrt(Lcont_sig)
    print 'NaID_alt Rref ', Rcont, np.sqrt(Rcont_sig)
    print 'NaID_alt ', NaID_alt, NaID_alt_sig
    print


    # HeI D3
    fig, (ax1, ax2, ax3) = plt.subplots(3)

    w_cent = 5875.620
    w_size = 0.2
    order =  53
    HeIcent, HeIcent_sig = line_flux(w_cent,w_size)
    #HeIcent, HeIcent_sig = line_integ(w_cent,w_size)
    add_subplot(ax1,'HeIcent',True)

    ## Reference region

    w_cent = 5869.000
    w_size = 2.5
    order =  53
    Lcont, Lcont_sig = line_flux(w_cent,w_size)
    #Lcont, Lcont_sig = line_integ(w_cent,w_size)
    add_subplot(ax2,'Left cont.')

    w_cent = 5881.000
    w_size = 2.5
    order =  53
    Rcont, Rcont_sig = line_flux(w_cent,w_size)
    #Rcont, Rcont_sig = line_integ(w_cent,w_size)
    add_subplot(ax3,'Right cont.')

    HeI = 2.* (HeIcent)/(Lcont+Rcont)
    HeI_sig = HeI * np.sqrt(HeIcent_sig/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)

    print 'HeI cent ', HeIcent, np.sqrt(HeIcent_sig)
    print 'HeI Lref ', Lcont, np.sqrt(Lcont_sig)
    print 'HeI Rref ', Rcont, np.sqrt(Rcont_sig)
    print 'HeI ', HeI, HeI_sig
    print
    fig.savefig(plot_dir + '/HeI_' + input_file + '.png', dpi = 300)


    fileout.write('{0:8s}\t{1:8f}\t{2:8f}\t{3:8f}\t{4:8f}\t{5:8f}\t{6:8f} '
    '{7:8f}\t{8:8f}\t{9:8f}\t{10:8f}\t{11:8f}\t{12:8f}\t{13:8f}\t{14:8f}\t{15:8f}\n'.format(\
       input_file,bjd,rv,rvn,R_out,R_sig,Halpha,Halpha_sig,
       NaID, NaID_sig, HeI, HeI_sig,
       HalphaHF, HalphaHF_sig, CaIquiet, CaIquiet_sig))


    s1d_hdu.close()

    fileout.close()
    #plt.show()
