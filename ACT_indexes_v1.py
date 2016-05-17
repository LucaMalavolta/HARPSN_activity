from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


def line_flux(w_cent,w_size,order):
   cond = (np.abs(wave_rest[order,:]-w_cent)<w_size)
   denom = sum(dwave_rest[order,cond])
   cont = sum(data_cor[order,cond])/denom
   cont_sig  = sum((data[order,cond]+noise**2)/(blaze[order,cond])**2)/denom**2
   return cont, cont_sig

def line_flux_triangle(w_cent,w_size,order):
   cond = (np.abs(wave_rest[order,:]-w_cent)<w_size)
   triangle = 1. - np.abs(wave_rest[order,:]-w_cent)/w_size
   denom = sum(triangle[cond]*dwave_rest[order,cond])
   cent = sum((data_cor[order,cond])*triangle[cond])/denom
   cent_sig = sum(  (data[order,cond]+noise**2)*(triangle[cond]/blaze[order,cond])**2 ) /denom**2
   return cent, cent_sig

print 'DISCLAIMER: the value obtained for the NaID are not good for star-to-star '
print 'comparison, since the continuum regions include several spectral lines '
print 'A different approach than Gomes da Silva 2011+ was followed'
print 'These value are meant only to assess the internal variation of the activity level'
print 'of a star'


#archive_dir = '/Users/malavolta/Astro/CODE/PyIndex'
archive_dir = '/media/malavolta/data/HARPN/data/'
input_list=np.genfromtxt('input.list',dtype='string')

ii_list = 0

cc = 299792.458

fileout = open('output.list','w')
fileout.write('#name 0, bjd 1 ,rv 2,rvn 3,R_out 4,R_sig 5,Halpha 6,Halpha_sig 7,HalphaHF 8,HalphaHF_sig 9,CaIquiet 10,CaIquiet_sig 11,NaID 12, NaID_sig 13,HeI 14,HeI_sig 15 \n')

#for ii_list in xrange(0,10):

for ii_list in xrange(0,np.size(input_list[:,0])):
   e2ds_file = archive_dir + '/' + input_list[ii_list,0] + '/' + input_list[ii_list,1] + '_e2ds_A.fits'


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
   blaze_file = archive_dir + '/' + input_list[ii_list,0] + '/' + e2ds_header['HIERARCH TNG DRS BLAZE FILE']
   blaze_hdu = fits.open(blaze_file)
   blaze = blaze_hdu[0].data
   data_cor = data/blaze

   ccf_file = archive_dir + '/' + input_list[ii_list,0] + '/' + input_list[ii_list,1] + '_ccf_' + input_list[ii_list,2] + '_A.fits'
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


   ## R H&K calculation

   w_cent = 3933.664
   w_size = 1.0900
   order =  1
   CaIIK, CaIIK_sig = line_flux_triangle(w_cent,w_size,order)
   print CaIIK, CaIIK_sig

   w_cent = 3968.470
   w_size = 1.0900
   order =  3
   CaIIH, CaIIH_sig = line_flux_triangle(w_cent,w_size,order)

   w_cent = 3901.070
   w_size = 10.
   order =  0
   Vcont, Vcont_sig = line_flux(w_cent,w_size,order)

   w_cent = 4001.070
   w_size = 10.
   order =  4
   Rcont, Rcont_sig = line_flux(w_cent,w_size,order)

   R_raw=(CaIIK+CaIIH)/(Vcont+Rcont)
   R_out = 1.111 * R_raw + 0.0153
   R_sig = 1.111 * np.sqrt((CaIIK_sig+CaIIH_sig)/(Vcont+Rcont)**2 + (CaIIK+CaIIH)**2*(Vcont_sig+Rcont_sig)/(Vcont+Rcont)**4)

   print 'CaII K ', CaIIK, np.sqrt(CaIIK_sig)
   print 'CaII H ', CaIIH, np.sqrt(CaIIH_sig)
   print 'CaII V ', Vcont, np.sqrt(Vcont_sig)
   print 'CaII R ', Rcont, np.sqrt(Rcont_sig)
   print 'R out  ', R_out, R_sig
   print


   ## Halpha Central
   ## Central wavelength at 6562.828 (Gomez da Silva 2011+, Robertson 2013a+)
   w_cent = 6562.828
   w_size = 0.8
   order =  64
   Hcent, Hcent_sig = line_flux(w_cent,w_size,order)

   # first continuum window in the range [-700:-300] km/s, as in Robertson 2013a+, Kurster
   w_cent = 6551.88238
   w_size = 4.37825
   order =  64
   Lcont, Lcont_sig = line_flux(w_cent,w_size,order)

   w_cent = 6580.341
   w_size = 4.37825
   order =  64
   Rcont, Rcont_sig = line_flux(w_cent,w_size,order)

   Halpha =(Hcent)/(Lcont+Rcont)
   Halpha_sig = Halpha * np.sqrt(Hcent_sig/Hcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
   print input_list[ii_list,1], input_list[ii_list,4]
   print 'H cent ', Hcent, np.sqrt(Hcent_sig)
   print 'H Lref ', Lcont, np.sqrt(Lcont_sig)
   print 'H Rref ', Rcont, np.sqrt(Rcont_sig)
   print 'Halpha ', Halpha, Halpha_sig
   print

   ## Halpha Central with half-size window
   ## Central wavelength at 6562.828 (Gomez da Silva 2011+, Robertson 2013a+)
   w_cent = 6562.828
   w_size = 0.4
   order =  64
   Hcent, Hcent_sig = line_flux(w_cent,w_size,order)

   HalphaHF =(Hcent)/(Lcont+Rcont)
   HalphaHF_sig = HalphaHF * np.sqrt(Hcent_sig/Hcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
   print input_list[ii_list,1], input_list[ii_list,4]
   print 'HalphaHF ', HalphaHF, HalphaHF_sig
   print




   ############### CaI activity-less index
   ## CaI "quiet" index (Kurster 2003+, Robertson 2013a+)
   w_cent = 6572.795
   w_size = 0.34
   order =  64
   CaIcent, CaIcent_sig = line_flux_triangle(w_cent,w_size,order)

   ## Rcont, Lcont as in the previous definition
   CaIquiet =(CaIcent)/(Lcont+Rcont)
   CaIquiet_sig = CaIquiet * np.sqrt(CaIcent_sig/CaIcent**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
   print 'CaI quiet ', CaIquiet, CaIquiet_sig
   print


   ############### NaI D
   ## activity index based on the flux in the cores of the Na I D1 & D2 lines
   ## as described in Sect. 3 of Gomes da Silva 2011+ and Dias 2007a+
   ## DISCLAIMER: the value obtained below is NOT good

   # NaI D1
   w_cent = 5895.920
   w_size = 0.25
   order =  53
   NaID1, NaID1_sig = line_flux(w_cent,w_size,order)

   w_cent = 5889.950
   w_size = 0.25
   order =  53
   NaID2, NaID2_sig = line_flux(w_cent,w_size,order)

   ## Reference region

   w_cent = 5805.000
   w_size = 5.0
   order =  52
   Lcont, Lcont_sig = line_flux(w_cent,w_size,order)

   w_cent = 6090.000
   w_size = 10.00
   order =  57
   Rcont, Rcont_sig = line_flux(w_cent,w_size,order)

   NaID =(NaID1+NaID2)/(Lcont+Rcont)
   NaID_sig = NaID * np.sqrt((NaID1_sig+NaID2_sig)/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
   print input_list[ii_list,1], input_list[ii_list,4]
   print 'NaID1 cen ', NaID1, np.sqrt(NaID1_sig)
   print 'NaID2 cen ', NaID2, np.sqrt(NaID2_sig)
   print 'NaID Lref ', Lcont, np.sqrt(Lcont_sig)
   print 'NaID Rref ', Rcont, np.sqrt(Rcont_sig)
   print 'NaID ', NaID, NaID_sig
   print


   # HeI D1
   w_cent = 5875.620
   w_size = 0.2
   order =  53
   HeIcent, HeIcent_sig = line_flux(w_cent,w_size,order)

   ## Reference region

   w_cent = 5869.000
   w_size = 2.5
   order =  53
   Lcont, Lcont_sig = line_flux(w_cent,w_size,order)

   w_cent = 5881.000
   w_size = 2.50
   order =  53
   Rcont, Rcont_sig = line_flux(w_cent,w_size,order)

   HeI =(HeIcent)/(Lcont+Rcont)
   HeI_sig = HeI * np.sqrt(HeIcent_sig/(NaID1+NaID2)**2 + (Lcont_sig+Rcont_sig)/(Lcont+Rcont)**2)
   print input_list[ii_list,1], input_list[ii_list,4]
   print 'HeI cent ', HeIcent, np.sqrt(HeIcent_sig)
   print 'HeI Lref ', Lcont, np.sqrt(Lcont_sig)
   print 'HeI Rref ', Rcont, np.sqrt(Rcont_sig)
   print 'HeI ', HeI, HeI_sig
   print



   #if ii_list==0:
   #  plt.xlim(5800,5900)
   #  for ii in xrange(51,66):
   #    plt.plot(wave_rest[ii,:],data_cor[ii,:])

   fileout.write('{0:8s} {1:8f} {2:8f} {3:8f} {4:8f} {5:8f} {6:8f} {7:8f} {8:8f} {9:8f} {10:8f} {11:8f} {12:8f} {13:8f} {14:8f} {15:8f} \n'.format(\
      input_list[ii_list,1],bjd,rv,rvn,R_out,R_sig,Halpha,Halpha_sig,HalphaHF,HalphaHF_sig,CaIquiet,CaIquiet_sig,NaID, NaID_sig,HeI,HeI_sig))


   e2ds_hdu.close()
   ccf_hdu.close()
   blaze_hdu.close()
fileout.close()
#plt.show()
