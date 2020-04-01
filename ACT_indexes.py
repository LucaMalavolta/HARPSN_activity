#! /usr/bin/python
# -*- coding: utf-8 -*-
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
#import seaborn as sns

__version__ = '2020-03-25'
__author__ = 'Luca Malavolta'

def get_instrument_keywords(obs_name):

    HARPN_dictionary = {
               'intrument_name' : 'HARPS-N',
               'deg_ll'  :   'HIERARCH TNG DRS CAL TH DEG LL',
               'naxis1' : 'NAXIS1',
               'naxis2' : 'NAXIS2',
               'sigdet' : 'HIERARCH TNG DRS CCD SIGDET',
               'gain'   : 'HIERARCH TNG DRS CCD CONAD',
               'coeff_ll': 'HIERARCH TNG DRS CAL TH COEFF LL',
               'blaze': 'HIERARCH TNG DRS BLAZE FILE',

                'bjd' : 'HIERARCH TNG DRS BJD',
                'berv': 'HIERARCH TNG DRS BERV',
                'rv'  : 'HIERARCH TNG DRS CCF RVC',
                'rvn' : 'HIERARCH TNG DRS CCF NOISE',

                'CaII_K': [3933.664, 1.0900, 1], #order in Python standard
                'CaII_H': [3968.470, 1.0900, 3],
                'CaII_V': [3901.070, 10.000, 0],
                'CaII_R': [4001.070, 10.000, 4],

                ## Halpha Central
                ## Central wavelength at 6562.828 (Gomez da Silva 2011+, Robertson 2013a+)
                'Halpha_cent': [6562.808, 0.8000, 64],
                # first continuum window in the range [-700:-300] km/s, as in Robertson 2013a+, Kurster
                'Halpha_Lcont': [6550.887, 5.3750, 64],
                'Halpha_Rcont': [6580.310, 4.3782, 64],

                ## Halpha Central with half-size window
                ## Central wavelength at 6562.828 (Gomez da Silva 2011+, Robertson 2013a+)
                'Halpha_half': [6562.808, 0.4000, 64],

                ############### CaI activity-less index
                ## CaI "quiet" index (Kurster 2003+, Robertson 2013a+)
                'CaI_cent': [6572.795, 0.3400, 64], #order in Python standard
                # first continuum window in the range [-700:-300] km/s, as in Robertson 2013a+, Kurster
                'CaI_Lcont': [6550.887, 5.3750, 64],
                'CaI_Rcont': [6580.310, 4.3782, 64],


                ############### NaI D alt
                ## Alternative approach to obtain Na I D1 & D2 lines
                ## check the disclaimer
                'NaI_D1': [5895.920, 0.2500, 53], #order in Python standard
                'NaI_D2': [5889.950, 0.2500, 53],
                'NaI_Lcont': [5805.000, 5.0000, 52],
                'NaI_Rcont': [6090.000, 10.000, 57],

                ############### NaI D
                ## NaID computation
                ## as described in Sect. 3 of Gomes da Silva 2011+ and Dias 2007a+
                ## check the disclaimer
                ## it uses the highest 10 instead on the flux


                'HeI_cent'   : [5875.620, 0.200, 53], #order in Python standard
                'HeI_Lcont'  : [5869.000, 2.500, 53],
                'HeI_Rcont'  : [5881.000, 2.500, 53],

    }

    HARPS_dictionary = {
               'intrument_name' : 'HARPS',
               'deg_ll'  :   'HIERARCH ESO DRS CAL TH DEG LL',
               'naxis1' : 'NAXIS1',
               'naxis2' : 'NAXIS2',
               'sigdet' : 'HIERARCH ESO DRS CCD SIGDET',
               'gain'   : 'HIERARCH ESO DRS CCD CONAD',
               'coeff_ll': 'HIERARCH ESO DRS CAL TH COEFF LL',
               'blaze': 'HIERARCH ESO DRS BLAZE FILE',

                'bjd' : 'HIERARCH ESO DRS BJD',
                'berv': 'HIERARCH ESO DRS BERV',
                'rv'  : 'HIERARCH ESO DRS CCF RVC',
                'rvn' : 'HIERARCH ESO DRS CCF NOISE',

                'CaII_K': [3933.664, 1.0900, 5], #order in Python standard
                'CaII_H': [3968.470, 1.0900, 7],
                'CaII_V': [3901.070, 10.000, 4],
                'CaII_R': [4001.070, 10.000, 8],

                ## Halpha Central
                ## Central wavelength at 6562.828 (Gomez da Silva 2011+, Robertson 2013a+)
                'Halpha_cent': [6562.808, 0.8000, 67],
                # first continuum window in the range [-700:-300] km/s, as in Robertson 2013a+, Kurster
                'Halpha_Lcont': [6550.887, 5.3750, 67],
                'Halpha_Rcont': [6580.310, 4.3782, 67],

                ## Halpha Central with half-size window
                ## Central wavelength at 6562.828 (Gomez da Silva 2011+, Robertson 2013a+)
                'Halpha_half': [6562.808, 0.4000, 67],

                ############### CaI activity-less index
                ## CaI "quiet" index (Kurster 2003+, Robertson 2013a+)
                'CaI_cent': [6572.795, 0.3400, 67], #order in Python standard
                # first continuum window in the range [-700:-300] km/s, as in Robertson 2013a+, Kurster
                'CaI_Lcont': [6550.887, 5.3750, 67],
                'CaI_Rcont': [6580.310, 4.3782, 67],


                ############### NaI D alt
                ## Alternative approach to obtain Na I D1 & D2 lines
                ## check the disclaimer
                'NaI_D1': [5895.920, 0.2500, 56], #order in Python standard
                'NaI_D2': [5889.950, 0.2500, 56],
                'NaI_Lcont': [5805.000, 5.0000, 55],
                'NaI_Rcont': [6090.000, 10.000, 59],  #vicino al bordo

                ############### NaI D
                ## NaID computation
                ## as described in Sect. 3 of Gomes da Silva 2011+ and Dias 2007a+
                ## check the disclaimer
                ## it uses the highest 10 instead on the flux


                'HeI_cent'   : [5875.620, 0.200, 56], #order in Python standard
                'HeI_Lcont'  : [5869.000, 2.500, 56],
                'HeI_Rcont'  : [5881.000, 2.500, 56],

    }


    print(obs_name)
    if 'HARPN' in obs_name:
        return HARPN_dictionary
    if 'HARPS' in obs_name:
        return HARPS_dictionary

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

def add_subplot(w_cent, w_size, order, plot_obj, title, wide='False', half_size=False):
   plot_obj.set_ylabel(title,fontsize=10)
   #plot_obj.title.set_text(title)

   # We change the fontsize of minor ticks label
   plot_obj.tick_params(axis='both', which='major', labelsize=10)
   plot_obj.tick_params(axis='both', which='minor', labelsize=8)

   if wide==True:
       print(title, np.amin(wave_rest[order,:]))
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

    print('DISCLAIMER: the value obtained for the NaID are not good for star-to-star comparison,')
    print('since the continuum regions include several spectral lines ')
    print('A different approach than Gomes da Silva 2011+ was followed')
    print('These value are meant only to assess the internal variation of the activity level of a star')


    parser = argparse.ArgumentParser(description='Generalized Lomb-Scargle periodogram.', add_help=False)
    argadd = parser.add_argument   # function short cut
    argadd('-?', '-h', '-help', '--help', help='show this help message and exit', action='help')
    argadd('file_list', nargs='?', help='Data file. If not specified an example will be shown.')
    argadd('archive_dir', nargs='?',default='./archive', help='Archive directory, default is "./archive" .')
    argadd('-nodate', nargs='?', default=False, help='Save full corellation plot - it may be slow!')
    args = parser.parse_args().__dict__
    file_list = args.pop('file_list')
    archive_dir = args.pop('archive_dir')
    nodate = args.pop('nodate')

    print(nodate)

    if file_list is None:
    # No data file given. Show example:
        print()
        print('Input file must have this structure:')
        print('2016-01-07   HARPN.2016-01-08T02-30-21.236   G2   A   TARG')
        print()
        quit()

    file_rad, file_ext = os.path.splitext(file_list)

    input_list=np.genfromtxt(file_list, dtype='U')

    plot_dir = './' + file_rad + '_plots'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)


    cc = 299792.458

    fileout = open(file_rad+'.output','w')
    fileout.write('file_name                    \tbjd           \trv      \trv_noise\tR_out   \tR_sig   \tHalpha  \tHa_sig  \tNaID    \tNaID_sig\tHeI     \tHeI_sig \tHalphaHF\tHaHF_sig\tCaIquiet\tCaIq_sig\n')
    fileout.write('#############################\t##############\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\t########\n')
    #for ii_list in range(0,10):

    #for input_night in input_list:
    for input_night   in input_list:

       if nodate is  False:
           input_rad = input_night[0] + '/'
       else:
           input_rad = ''

       instrument_dict = get_instrument_keywords(input_night[1])
       print('Analysis of:', input_night)
       print('Instrument: ',instrument_dict['intrument_name'] )
       print()

       e2ds_file = archive_dir +  '/' + input_rad + input_night[1] + '_e2ds_A.fits'


       e2ds_hdu = fits.open(e2ds_file)
       e2ds_header = e2ds_hdu[0].header
       data = e2ds_hdu[0].data
       deg  = e2ds_header[instrument_dict['deg_ll']]
       naxis1 = e2ds_header[instrument_dict['naxis1']]
       naxis2 = e2ds_header[instrument_dict['naxis2']]

       sigdet = e2ds_header[instrument_dict['sigdet']]
       gain   = e2ds_header[instrument_dict['gain']]
       noise = sigdet*np.sqrt(6)*gain

       #x = np.arange(naxis1-1,-1,-1.)
       x = np.arange(0,naxis1)
       wave = np.zeros([naxis2,naxis1],dtype=np.double)

       for n in range(0,naxis2):
          for i in range(deg,-1,-1):
           a_sel = i + n*(1+deg)
           a_coeff = e2ds_header[instrument_dict['coeff_ll']+repr(a_sel)]
           if (i == deg):
           	wave[n,:] = a_coeff
           else:
           	wave[n,:] = wave[n,:]*x + a_coeff

       #Extracting keyword info
       blaze_file = archive_dir +  '/' + input_rad + e2ds_header[instrument_dict['blaze']]
       blaze_hdu = fits.open(blaze_file)
       blaze = blaze_hdu[0].data
       data_cor = data/blaze

       ccf_file = archive_dir +  '/' + input_rad + input_night[1] + '_ccf_' + input_night[2] + '_A.fits'
       ccf_hdu = fits.open(ccf_file)
       ccf_header = ccf_hdu[0].header


       print('e2ds file:   ', e2ds_file)
       print('ccf file:    ', ccf_file)
       print('blaze file:  ', blaze_file)
       print()

       bjd  = ccf_header[instrument_dict['bjd']]
       berv = ccf_header[instrument_dict['berv']]
       rv   = ccf_header[instrument_dict['rv']]
       rvn  = ccf_header[instrument_dict['rvn']]

       wave_rest = wave * ((1.+berv/cc)/(1.+rv/cc))
       dwave_rest = np.zeros([naxis2,naxis1],dtype=np.double)
       dwave_rest[:,1:] = wave_rest[:,1:]-wave_rest[:,:-1]
       dwave_rest[:,0] = dwave_rest[:,1]

       ## R H&K calculation
       fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)
       order_ref = instrument_dict['CaII_R'][-1]
       norm_factor = np.sum(data_cor[order_ref,:])/np.sum(dwave_rest[order_ref,:])
       #plots are normalized acording to the flux of the chosen order

       CaIIK, CaIIK_sig = line_flux_triangle(*instrument_dict['CaII_K'])
       add_subplot(*instrument_dict['CaII_K'], ax1,'CaII K',True)

       CaIIH, CaIIH_sig = line_flux_triangle(*instrument_dict['CaII_H'])
       add_subplot(*instrument_dict['CaII_H'], ax2,'CaII H',True)

       Vcont, Vcont_sig = line_flux(*instrument_dict['CaII_V'])
       add_subplot(*instrument_dict['CaII_V'], ax3,'Left cont.')

       Rcont, Rcont_sig = line_flux(*instrument_dict['CaII_R'])
       add_subplot(*instrument_dict['CaII_R'], ax4,'Right cont.')

       R_raw=(CaIIK+CaIIH)/(Vcont+Rcont)
       R_out = 1.111 * R_raw + 0.0153
       R_sig = 1.111 * np.sqrt((CaIIK_sig+CaIIH_sig)/(Vcont+Rcont)**2 + (CaIIK+CaIIH)**2*(Vcont_sig+Rcont_sig)/(Vcont+Rcont)**4)

       print('CaII K ', CaIIK, np.sqrt(CaIIK_sig))
       print('CaII H ', CaIIH, np.sqrt(CaIIH_sig))
       print('CaII V ', Vcont, np.sqrt(Vcont_sig))
       print('CaII R ', Rcont, np.sqrt(Rcont_sig))
       print('R out  ', R_out, R_sig)
       print()
       fig.savefig(plot_dir + '/CaII_' + input_night[1] + '.png', dpi = 300)


       fig, (ax1, ax2, ax3) = plt.subplots(3)

       order_ref = instrument_dict['Halpha_Rcont'][-1]
       norm_factor = np.sum(data_cor[order_ref,:])/np.sum(dwave_rest[order_ref,:])


       Hcent, Hcent_sig = line_flux(*instrument_dict['Halpha_cent'])
       add_subplot(*instrument_dict['Halpha_cent'],ax1,'Halpha cent',True,True)

       Lcont, Lcont_sig = line_flux(*instrument_dict['Halpha_Lcont'])
       add_subplot(*instrument_dict['Halpha_Lcont'],ax2,'Left cont.',False)

       Rcont, Rcont_sig = line_flux(*instrument_dict['Halpha_Rcont'])
       add_subplot(*instrument_dict['Halpha_Rcont'],ax3,'Right continuum',False)

       Halpha =(Hcent)/(Lcont+Rcont)
       Halpha_sig = Halpha * np.sqrt(Hcent_sig/Hcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print(input_night[1], input_night[4])
       print('H cent ', Hcent, np.sqrt(Hcent_sig))
       print('H Lref ', Lcont, np.sqrt(Lcont_sig))
       print('H Rref ', Rcont, np.sqrt(Rcont_sig))
       print('Halpha ', Halpha, Halpha_sig)
       print()

       HcentHF, HcentHF_sig = line_flux(*instrument_dict['Halpha_half'])
       add_subplot(*instrument_dict['Halpha_half'],ax1,'Halpha (half width)',True,True)

       HalphaHF =(HcentHF)/(Lcont+Rcont)
       HalphaHF_sig = HalphaHF * np.sqrt(HcentHF_sig/HcentHF**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print(input_night[1], input_night[4])
       print('HalphaHF ', HalphaHF, HalphaHF_sig)
       print()

       fig.savefig(plot_dir + '/Halpha_' + input_night[1] + '.png', dpi = 300)

       fig, (ax1, ax2, ax3) = plt.subplots(3)

       order_ref = instrument_dict['CaI_Rcont'][-1]
       norm_factor = np.sum(data_cor[order_ref,:])/np.sum(dwave_rest[order_ref,:])

       CaIcent, CaIcent_sig = line_flux_triangle(*instrument_dict['CaI_cent'])
       add_subplot(*instrument_dict['CaI_cent'],ax1,'CaI central',True,True)

       Lcont, Lcont_sig = line_flux(*instrument_dict['CaI_Lcont'])
       add_subplot(*instrument_dict['CaI_Lcont'],ax2,'Left continuum',False)

       Rcont, Rcont_sig = line_flux(*instrument_dict['CaI_Rcont'])
       add_subplot(*instrument_dict['CaI_Rcont'],ax3,'Right continuum',False)

       ## Rcont, Lcont as in the previous definition
       CaIquiet = 2.*(CaIcent)/(Lcont+Rcont)
       CaIquiet_sig = CaIquiet * np.sqrt(CaIcent_sig/CaIcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print('CaI quiet ', CaIquiet, CaIquiet_sig)
       print()
       fig.savefig(plot_dir + '/CaIquiet_' + input_night[1] + '.png', dpi = 300)




       fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)
       order_ref = instrument_dict['NaI_Rcont'][-1]
       norm_factor = np.sum(data_cor[order_ref,:])/np.sum(dwave_rest[order_ref,:])


       NaID1, NaID1_sig = line_flux(*instrument_dict['NaI_D1'])
       #NaID1, NaID1_sig = line_integ(w_cent,w_size,order)
       add_subplot(*instrument_dict['NaI_D1'],ax1,'NaI D1',True)

       NaID2, NaID2_sig = line_flux(*instrument_dict['NaI_D2'])
       #NaID2, NaID2_sig = line_integ(w_cent,w_size,order)
       add_subplot(*instrument_dict['NaI_D2'],ax2,'NaI D2',True)

       ## Reference region
       Lcont, Lcont_sig = line_flux(*instrument_dict['NaI_Lcont'])
       #Lcont, Lcont_sig = line_integ(*instrument_dict['NaI_Lcont'])
       add_subplot(*instrument_dict['NaI_Lcont'],ax3,'Left continuum')

       Rcont, Rcont_sig = line_flux(*instrument_dict['NaI_Rcont'])
       #Rcont, Rcont_sig = line_integ(*instrument_dict['NaI_Rcont'])
       add_subplot(*instrument_dict['NaI_Rcont'], ax4,'Right continuum')

       NaID =(NaID1+NaID2)/(Lcont+Rcont)
       NaID_sig = NaID * np.sqrt((NaID1_sig+NaID2_sig)/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print(input_night[1], input_night[4])
       print('NaI D1 cen ', NaID1, np.sqrt(NaID1_sig))
       print('NaI D2 cen ', NaID2, np.sqrt(NaID2_sig))
       print('NaI D Lref ', Lcont, np.sqrt(Lcont_sig))
       print('NaI D Rref ', Rcont, np.sqrt(Rcont_sig))
       print('NaI D ', NaID, NaID_sig)
       print()
       fig.savefig(plot_dir + '/NaID_' + input_night[1] + '.png', dpi = 300)

       ############### NaI D
       ## NaID computation
       ## as described in Sect. 3 of Gomes da Silva 2011+ and Dias 2007a+
       ## check the disclaimer

       NaID1, NaID1_sig = highest10(*instrument_dict['NaI_D1'])

       NaID2, NaID2_sig = highest10(*instrument_dict['NaI_D2'])

       ## Reference region
       Lcont, Lcont_sig = highest10(*instrument_dict['NaI_Lcont'])

       Rcont, Rcont_sig = highest10(*instrument_dict['NaI_Rcont'])


       NaID_alt =(NaID1+NaID2)/(Lcont+Rcont)
       NaID_alt_sig = NaID_alt * np.sqrt((NaID1_sig+NaID2_sig)/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print(input_night[1], input_night[4])
       print('NaID1_alt cen ', NaID1, np.sqrt(NaID1_sig))
       print('NaID2_alt cen ', NaID2, np.sqrt(NaID2_sig))
       print('NaID_alt Lref ', Lcont, np.sqrt(Lcont_sig))
       print('NaID_alt Rref ', Rcont, np.sqrt(Rcont_sig))
       print('NaID_alt ', NaID_alt, NaID_alt_sig)
       print()


       # HeI D3
       fig, (ax1, ax2, ax3) = plt.subplots(3)
       order_ref = instrument_dict['HeI_cent'][-1]
       norm_factor = np.sum(data_cor[order_ref,:])/np.sum(dwave_rest[order_ref,:])

       HeIcent, HeIcent_sig = line_flux(*instrument_dict['HeI_cent'])
       #HeIcent, HeIcent_sig = line_integ(*instrument_dict['HeI_cent'])
       add_subplot(*instrument_dict['HeI_cent'],ax1,'HeI cent',True)

       ## Reference region

       Lcont, Lcont_sig = line_flux(*instrument_dict['HeI_Lcont'])
       #Lcont, Lcont_sig = line_integ(*instrument_dict['HeI_Lcont'])
       add_subplot(*instrument_dict['HeI_Lcont'],ax2,'Left continumm')

       Rcont, Rcont_sig = line_flux(*instrument_dict['HeI_Rcont'])
       #Rcont, Rcont_sig = line_integ(*instrument_dict['HeI_Rcont'])
       add_subplot(*instrument_dict['HeI_Rcont'],ax3,'Right continuum')

       HeI = 2.* (HeIcent)/(Lcont+Rcont)
       HeI_sig = HeI * np.sqrt(HeIcent_sig/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
       print(input_night[1], input_night[4])
       print('HeI cent ', HeIcent, np.sqrt(HeIcent_sig))
       print('HeI Lref ', Lcont, np.sqrt(Lcont_sig))
       print('HeI Rref ', Rcont, np.sqrt(Rcont_sig))
       print('HeI ', HeI, HeI_sig)
       print()
       fig.savefig(plot_dir + '/HeI_' + input_night[1] + '.png', dpi = 300)

       fileout.write('{0:8s}\t{1:8f}\t{2:8f}\t{3:8f}\t{4:8f}\t{5:8f}\t{6:8f} '
       '{7:8f}\t{8:8f}\t{9:8f}\t{10:8f}\t{11:8f}\t{12:8f}\t{13:8f}\t{14:8f}\t{15:8f}\n'.format(\
          input_night[1],bjd,rv,rvn,R_out,R_sig,Halpha,Halpha_sig,
          NaID, NaID_sig, HeI, HeI_sig,
          HalphaHF, HalphaHF_sig, CaIquiet, CaIquiet_sig))


       e2ds_hdu.close()
       ccf_hdu.close()
       blaze_hdu.close()

       print(input_night, 'analysis completed')
       print()
       print()
    fileout.close()
    #plt.show()
