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

def line_flux(w_cent,w_size,order):
   cond = (np.abs(wave_rest[order,:]-w_cent)<=w_size)
   wl_range = sum(dwave_rest[order,cond])
   cont = sum(data_cor[order,cond]) / wl_range
   cont_sig  = sum((data[order,cond]+noise**2)/(blaze[order,cond])**2)/ wl_range**2
   return cont, cont_sig

def line_integ(w_cent,w_size,order):
   cond = (np.abs(wave_rest[order,:]-w_cent)<=w_size)
   wl_range = sum(dwave_rest[order,cond])
   cont = sum(data_cor[order,cond]) / wl_range * w_size * 2.
   cont_sig  = sum((data[order,cond]+noise**2)/(blaze[order,cond])**2) #* (w_size * 2. / wl_range)**2
   return cont, cont_sig

def highest10(w_cent,w_size,order):
   cond = (np.abs(wave_rest[order,:]-w_cent)<=w_size)
   vals = data_cor[order,cond]/dwave_rest[order,cond]
   ind = np.argpartition(vals, -10)[-10:]
   cont = np.average(vals[ind])
   cont_sig = (np.std(vals[ind])**2)/9.
   return cont, cont_sig

def line_flux_triangle(w_cent,w_size,order):
   cond = (np.abs(wave_rest[order,:]-w_cent)<=w_size)
   triangle = 1. - np.abs(wave_rest[order,:]-w_cent)/w_size
   denom = sum(triangle[cond]*dwave_rest[order,cond])
   cent = sum((data_cor[order,cond])*triangle[cond])/denom
   cent_sig = sum(  (data[order,cond]+noise**2)*(triangle[cond]/blaze[order,cond])**2 ) /denom**2
   return cent, cent_sig

def add_subplot(plot_obj, title, wide='False', half_size=False):
   plot_obj.set_ylabel(title,fontsize=10)
   #plot_obj.title.set_text(title)

   # We change the fontsize of minor ticks label
   plot_obj.tick_params(axis='both', which='major', labelsize=10)
   plot_obj.tick_params(axis='both', which='minor', labelsize=8)

   if wide==True:
       print title, np.amin(wave_rest[order,:])
       xlim_min = np.amax([w_cent-w_size*5., np.amin(wave_rest[order,:])])
       xlim_max = np.amin([w_cent+w_size*5., np.amax(wave_rest[order,:])])
   else:
       xlim_min = w_cent-w_size*2.
       xlim_max = w_cent+w_size*2.
   sel = (wave_rest[order,:]>xlim_min) & (wave_rest[order,:]<xlim_max)
   #ylim = np.max(data_cor[order,sel]/norm_factor/dwave_rest[order,sel])
   plot_obj.set_xlim(xlim_min,xlim_max)
   plot_obj.set_ylim(0.,1.5)
   plot_obj.axvline(w_cent-w_size/2.0, c='r')
   plot_obj.axvline(w_cent+w_size/2.0, c='r')
   if half_size:
      plot_obj.axvline(w_cent-w_size/4.0, c='g')
      plot_obj.axvline(w_cent+w_size/4.0, c='g')
   plot_obj.plot(wave_rest[order,sel],data_cor[order,sel]/norm_factor/dwave_rest[order,sel])

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
    argadd('archive_dir', nargs='?',help='Data file. If not specified example will be shown.')

    args = parser.parse_args().__dict__
    input_file = args.pop('input_file')
    archive_dir = args.pop('archive_dir')
    if input_file is None:
    # No data file given. Show example:
        print '2016-01-07   HARPN.2016-01-08T02-30-21.236   G2   A   TARG'
    if archive_dir is None:
    # No data file given. Show example:
        print 'using default archove directory ./archive'
        archive_dir = './archive'

    file_rad, file_ext = os.path.splitext(input_file)

    input_list=np.genfromtxt(input_file,dtype='string')

    plot_dir = './' + file_rad + '_plots'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)


    cc = 299792.458

    fileout = open(file_rad+'.output','w')
    fileout.write('file_name                    \tbjd           \trv      \trv_noise\tR_out   \tR_sig   \tHalpha  \tHa_sig  \tNaID    \tNaID_sig\tHeI     \tHeI_sig \tHalphaHF\tHaHF_sig\tCaIquiet\tCaIq_sig\n')
    fileout.write('#############################\t##############\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\n')
    #for ii_list in xrange(0,10):

    for input_night in input_list:

       e2ds_file = archive_dir +  '/' + input_night[1] + '_e2ds_A.fits'


       e2ds_hdu = fits.open(e2ds_file)
       e2ds_header = e2ds_hdu[0].header
       data = e2ds_hdu[0].data
       deg  = e2ds_header['HIERARCH TNG DRS CAL TH DEG LL']
       naxis1 = e2ds_header['NAXIS1']
       naxis2 = e2ds_header['NAXIS2']

       sigdet = e2ds_header['HIERARCH TNG DRS CCD SIGDET']
       gain   = e2ds_header['HIERARCH TNG DRS CCD CONAD']
       noise = sigdet*np.sqrt(6)*gain

       #x = np.arange(naxis1-1,-1,-1.)
       x = np.arange(0,naxis1)
       wave = np.zeros([naxis2,naxis1],dtype=np.double)

       for n in xrange(0,naxis2):
          for i in xrange(deg,-1,-1):
           a_sel = i + n*(1+deg)
           a_coeff = e2ds_header['HIERARCH TNG DRS CAL TH COEFF LL'+`a_sel`]
           if (i == deg):
           	wave[n,:] = a_coeff
           else:
           	wave[n,:] = wave[n,:]*x + a_coeff

       #Extracting keyword info
       blaze_file = archive_dir +  '/' + e2ds_header['HIERARCH TNG DRS BLAZE FILE']
       blaze_hdu = fits.open(blaze_file)
       blaze = blaze_hdu[0].data
       data_cor = data/blaze

       ccf_file = archive_dir +  '/' + input_night[1] + '_ccf_' + input_night[2] + '_A.fits'
       ccf_hdu = fits.open(ccf_file)
       ccf_header = ccf_hdu[0].header

       print e2ds_file
       print ccf_file
       print blaze_file

       bjd  = ccf_header['HIERARCH TNG DRS BJD']
       berv = ccf_header['HIERARCH TNG DRS BERV']
       rv   = ccf_header['HIERARCH TNG DRS CCF RVC']
       rvn  = ccf_header['HIERARCH TNG DRS CCF NOISE']

       wave_rest = wave * ((1.+berv/cc)/(1.+rv/cc))
       dwave_rest = np.zeros([naxis2,naxis1],dtype=np.double)
       dwave_rest[:,1:] = wave_rest[:,1:]-wave_rest[:,:-1]
       dwave_rest[:,0] = dwave_rest[:,1]

       ## R H&K calculation
       fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)
       order_ref = 4
       norm_factor = np.sum(data_cor[order_ref,:])/np.sum(dwave_rest[order_ref,:])
       #plots are normalized acording to the flux of the chosen order

       w_cent = 3933.664
       w_size = 1.0900
       order =  1
       CaIIK, CaIIK_sig = line_flux_triangle(w_cent,w_size,order)
       add_subplot(ax1,'CaII K',True)

       w_cent = 3968.470
       w_size = 1.0900
       order =  3
       CaIIH, CaIIH_sig = line_flux_triangle(w_cent,w_size,order)
       add_subplot(ax2,'CaII H',True)

       w_cent = 3901.070
       w_size = 10.
       order =  0
       Vcont, Vcont_sig = line_flux(w_cent,w_size,order)
       add_subplot(ax3,'Left cont.')

       w_cent = 4001.070
       w_size = 10.
       order =  4
       Rcont, Rcont_sig = line_flux(w_cent,w_size,order)
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
       fig.savefig(plot_dir + '/CaII_' + input_night[1] + '.png', dpi = 300)


       fig, (ax1, ax2, ax3) = plt.subplots(3)
       fig2, (ax4, ax5, ax6) = plt.subplots(3)

       order_ref = 64
       norm_factor = np.sum(data_cor[order_ref,:])/np.sum(dwave_rest[order_ref,:])

       ## Halpha Central
       ## Central wavelength at 6562.828 (Gomez da Silva 2011+, Robertson 2013a+)
       w_cent = 6562.808
       w_size = 0.8
       order =  64
       Hcent, Hcent_sig = line_flux(w_cent,w_size,order)
       add_subplot(ax1,'Halpha cent',True,True)

       # first continuum window in the range [-700:-300] km/s, as in Robertson 2013a+, Kurster
       w_cent = 6550.887
       w_size = 5.375
       order =  64
       Lcont, Lcont_sig = line_flux(w_cent,w_size,order)
       add_subplot(ax2,'Left cont.',False)
       add_subplot(ax5,'Left cont.',False)

       w_cent = 6580.31
       w_size = 4.37825
       order =  64
       Rcont, Rcont_sig = line_flux(w_cent,w_size,order)
       add_subplot(ax3,'Right continuum',False)
       add_subplot(ax6,'Right continuum',False)

       Halpha =(Hcent)/(Lcont+Rcont)
       Halpha_sig = Halpha * np.sqrt(Hcent_sig/Hcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print input_night[1], input_night[4]
       print 'H cent ', Hcent, np.sqrt(Hcent_sig)
       print 'H Lref ', Lcont, np.sqrt(Lcont_sig)
       print 'H Rref ', Rcont, np.sqrt(Rcont_sig)
       print 'Halpha ', Halpha, Halpha_sig
       print
       fig.savefig(plot_dir + '/Halpha_' + input_night[1] + '.png', dpi = 300)

       ## Halpha Central with half-size window
       ## Central wavelength at 6562.828 (Gomez da Silva 2011+, Robertson 2013a+)
       w_cent = 6562.808
       w_size = 0.4
       order =  64
       Hcent, Hcent_sig = line_flux(w_cent,w_size,order)

       HalphaHF =(Hcent)/(Lcont+Rcont)
       HalphaHF_sig = HalphaHF * np.sqrt(Hcent_sig/Hcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print input_night[1], input_night[4]
       print 'HalphaHF ', HalphaHF, HalphaHF_sig
       print



       ############### CaI activity-less index
       ## CaI "quiet" index (Kurster 2003+, Robertson 2013a+)


       order_ref = 64
       norm_factor = np.sum(data_cor[order_ref,:])/np.sum(dwave_rest[order_ref,:])

       w_cent = 6572.795
       w_size = 0.34
       order =  64
       CaIcent, CaIcent_sig = line_flux_triangle(w_cent,w_size,order)

       ## Rcont, Lcont as in the previous definition
       CaIquiet = 2.*(CaIcent)/(Lcont+Rcont)
       CaIquiet_sig = CaIquiet * np.sqrt(CaIcent_sig/CaIcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print 'CaI quiet ', CaIquiet, CaIquiet_sig
       print
       add_subplot(ax4,'CaI central',True)
       fig2.savefig(plot_dir + '/CaIquiet_' + input_night[1] + '.png', dpi = 300)



       ############### NaI D alt
       ## Alternative approach to obtain Na I D1 & D2 lines
       ## check the disclaimer
       fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)
       order_ref = 53
       norm_factor = np.sum(data_cor[order_ref,:])/np.sum(dwave_rest[order_ref,:])

       w_cent = 5895.920
       w_size = 0.25
       order =  53
       NaID1, NaID1_sig = line_flux(w_cent,w_size,order)
       #NaID1, NaID1_sig = line_integ(w_cent,w_size,order)
       add_subplot(ax1,'NaID1',True)

       w_cent = 5889.950
       w_size = 0.25
       order =  53
       NaID2, NaID2_sig = line_flux(w_cent,w_size,order)
       #NaID2, NaID2_sig = line_integ(w_cent,w_size,order)
       add_subplot(ax2,'NaID2',True)

       ## Reference region
       w_cent = 5805.000
       w_size = 5.0
       order =  52
       Lcont, Lcont_sig = line_flux(w_cent,w_size,order)
       #Lcont, Lcont_sig = line_integ(w_cent,w_size,order)
       add_subplot(ax3,'Left cont.')

       w_cent = 6090.000
       w_size = 10.00
       order =  57
       Rcont, Rcont_sig = line_flux(w_cent,w_size,order)
       #Rcont, Rcont_sig = line_integ(w_cent,w_size,order)
       add_subplot(ax4,'Right cont.')


       NaID =(NaID1+NaID2)/(Lcont+Rcont)
       NaID_sig = NaID * np.sqrt((NaID1_sig+NaID2_sig)/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print input_night[1], input_night[4]
       print 'NaID1 cen ', NaID1, np.sqrt(NaID1_sig)
       print 'NaID2 cen ', NaID2, np.sqrt(NaID2_sig)
       print 'NaID Lref ', Lcont, np.sqrt(Lcont_sig)
       print 'NaID Rref ', Rcont, np.sqrt(Rcont_sig)
       print 'NaID ', NaID, NaID_sig
       print
       fig.savefig(plot_dir + '/NaID_' + input_night[1] + '.png', dpi = 300)

       ############### NaI D
       ## NaID computation
       ## as described in Sect. 3 of Gomes da Silva 2011+ and Dias 2007a+
       ## check the disclaimer

       w_cent = 5895.920
       w_size = 0.25
       order =  53
       NaID1, NaID1_sig = highest10(w_cent,w_size,order)

       w_cent = 5889.950
       w_size = 0.25
       order =  53
       NaID2, NaID2_sig = highest10(w_cent,w_size,order)

       ## Reference region

       w_cent = 5805.000
       w_size = 5.0
       order =  52
       Lcont, Lcont_sig = highest10(w_cent,w_size,order)

       w_cent = 6090.000
       w_size = 10.00
       order =  57
       Rcont, Rcont_sig = highest10(w_cent,w_size,order)


       NaID_alt =(NaID1+NaID2)/(Lcont+Rcont)
       NaID_alt_sig = NaID_alt * np.sqrt((NaID1_sig+NaID2_sig)/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print input_night[1], input_night[4]
       print 'NaID1_alt cen ', NaID1, np.sqrt(NaID1_sig)
       print 'NaID2_alt cen ', NaID2, np.sqrt(NaID2_sig)
       print 'NaID_alt Lref ', Lcont, np.sqrt(Lcont_sig)
       print 'NaID_alt Rref ', Rcont, np.sqrt(Rcont_sig)
       print 'NaID_alt ', NaID_alt, NaID_alt_sig
       print


       # HeI D3
       fig, (ax1, ax2, ax3) = plt.subplots(3)
       order_ref = 53
       norm_factor = np.sum(data_cor[order_ref,:])/np.sum(dwave_rest[order_ref,:])

       w_cent = 5875.620
       w_size = 0.2
       order =  53
       HeIcent, HeIcent_sig = line_flux(w_cent,w_size,order)
       #HeIcent, HeIcent_sig = line_integ(w_cent,w_size,order)
       add_subplot(ax1,'HeIcent',True)

       ## Reference region

       w_cent = 5869.000
       w_size = 2.5
       order =  53
       Lcont, Lcont_sig = line_flux(w_cent,w_size,order)
       #Lcont, Lcont_sig = line_integ(w_cent,w_size,order)
       add_subplot(ax2,'Left cont.')

       w_cent = 5881.000
       w_size = 2.5
       order =  53
       Rcont, Rcont_sig = line_flux(w_cent,w_size,order)
       #Rcont, Rcont_sig = line_integ(w_cent,w_size,order)
       add_subplot(ax3,'Right cont.')

       HeI = 2.* (HeIcent)/(Lcont+Rcont)
       HeI_sig = HeI * np.sqrt(HeIcent_sig/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print input_night[1], input_night[4]
       print 'HeI cent ', HeIcent, np.sqrt(HeIcent_sig)
       print 'HeI Lref ', Lcont, np.sqrt(Lcont_sig)
       print 'HeI Rref ', Rcont, np.sqrt(Rcont_sig)
       print 'HeI ', HeI, HeI_sig
       print
       fig.savefig(plot_dir + '/HeI_' + input_night[1] + '.png', dpi = 300)


       fileout.write('{0:8s}\t{1:8f}\t{2:8f}\t{3:8f}\t{4:8f}\t{5:8f}\t{6:8f} '
       '{7:8f}\t{8:8f}\t{9:8f}\t{10:8f}\t{11:8f}\t{12:8f}\t{13:8f}\t{14:8f}\t{15:8f}\n'.format(\
          input_night[1],bjd,rv,rvn,R_out,R_sig,Halpha,Halpha_sig,
          NaID, NaID_sig, HeI, HeI_sig,
          HalphaHF, HalphaHF_sig, CaIquiet, CaIquiet_sig))


       e2ds_hdu.close()
       ccf_hdu.close()
       blaze_hdu.close()
    fileout.close()
    #plt.show()
