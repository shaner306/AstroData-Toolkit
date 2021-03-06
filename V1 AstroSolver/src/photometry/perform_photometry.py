# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 15:05:05 2022

@author: mstew

This file stores code used for photometry general_tools

Breakdown of this module:

┌───────────┐        ┌────────────────────────────────┐    ┌────────────────┐
│  Inputs   │        │  Phot. Methods:                │    │ Phot Results   │
├───────────┤        │ ┌────────────┐                 │    ├────────────────┤
├───────────┤        │ │PSF Fitting │                 │    ├────────────────┤
│           │        │ └────────────┘                 │    │                │
│Source     │        │ OR                             │    │X,Y pix location│
│  Estimates│        │ ┌────────────┐                 │    │                │
│           │        │ │Aperture    │                 │    │flux            │
│FWHM       │        │ └─────┬──────┘                 │    │                │
│           │        │       │                        │    │flux error      │
│imgdata    ├───────►│       │       ┌─────────────┐  ├───►│  (1 sigma)     │
│           │        │       ├──────►│Mean Aperture│  │    │                │
│bkg        │        │       │       └─────────────┘  │    │                │
│           │        │       │ OR                     │    │                │
│filepath   │        │       │       ┌─────────────┐  │    │                │
│           │        │       └──────►│Dynamic      │  │    │                │
│fits. hdr  │        │               └─────────────┘  │    │                │
│           │        │                                │    │                │
└───────────┘        └────────────────────────────────┘    └────────────────┘



"""
import numpy as np
from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import SigmaClip
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.table import Table
from astropy.io import fits
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF
import numpy as np
from photutils.aperture import CircularAperture,aperture_photometry,ApertureStats,CircularAnnulus
from photutils.utils import calc_total_error


def perform_PSF_photometry(irafsources, fwhm, imgdata, bkg,filepath,hdr,
                       fitter=LevMarLSQFitter(), fitshape=5):
    """
    Perform PSF photometry on all sources in a selected image.

    Parameters
    ----------
    irafsources : astropy.table.Table
        Table containing information of all stars detected in the image.
        Has columns:
            id,
            xcentroid,
            ycentroid,
            fwhm,
            sharpness,
            roundness,
            pa,
            npix,
            sky,
            peak,
            flux,
            mag
    fwhm : float
        Mean FWHM of all sources in the image.
    imgdata : numpy.ndarray
        Data from the fits file.
    bkg : float, optional
        Background value of the image in ADU. The default is None.
    bkg_estimator : callable, instance of any
    photutils.background.BackgroundBase subclass, optional
        bkg_estimator should be able to compute either a scalar background or
        a 2D background of a given 2D image.
        If None, no background subtraction is performed. The default is None.
    fitter : astropy.modeling.fitting.Fitter instance, optional
        Fitter object used to compute the optimized centroid positions and/or
        flux of the identified sources.
        The default is LevMarLSQFitter().
    fitshape : int or length-2 array-like, optional
        Rectangular shape around the center of a star which will be used to
        collect the data to do the fitting.
        Can be an integer to be the same along both axes. For example, 5 is
        the same as (5, 5),
        which means to fit only at the following relative pixel positions:
            [-2, -1, 0, 1, 2].
        Each element of fitshape must be an odd number. The default is 25.

    Returns
    -------
    photometry_result : astropy.table.Table or None
        Table with the photometry results, i.e., centroids and fluxes
        estimations and the initial estimates used to
        start the fitting process. Uncertainties on the fitted parameters
        are reported as columns called
        <paramname>_unc provided that the fitter object contains a dictionary
        called fit_info with the key param_cov,
        which contains the covariance matrix. If param_cov is not present,
        uncertanties are not reported.

    """
    daogroup = DAOGroup(2 * fwhm)
    psf_model = IntegratedGaussianPRF(sigma=fwhm * gaussian_fwhm_to_sigma)
    psf_model.x_0.fixed = True
    psf_model.y_0.fixed = True
    pos = Table(names=['x_0', 'y_0', 'flux_0'],
                data=[irafsources['xcentroid'],
                      irafsources['ycentroid'],
                      irafsources['flux']])

    photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=None,
                                    psf_model=psf_model,
                                    fitter=fitter,
                                    fitshape=fitshape)
    
    photometry_result = photometry(image=imgdata - bkg, init_guesses=pos)
    
    
    #Get Residual Image
    # residual_image=photometry.get_residual_image()
    #
    # fits_reisdual_image=fits.PrimaryHDU(data=residual_image,header=hdr)
    # fits_reisdual_image.writeto((filepath.split('.fits')[0]+'_residual.fits'),overwrite=True)
    return photometry_result


def perform_PSF_photometry_2(irafsources, fwhm, imgdata, bkg,filepath,hdr,
                       fitter=LevMarLSQFitter(), fitshape=5,produce_residual_data=False):
    
    #TODO: Find Dynamic way to adjust the PSF fitting

    '''
    An experimental version of the 'perform_PSF_photometry' function which uses
    3*fwhm of the chosen source to create a fit shape of the object.
    The primary reason for doing this is to explore the consequences
    of changing the fit shape of the basic psf function
    
    '''
    pos = Table(names=['x_0', 'y_0', 'flux_0'],
                data=[irafsources['xcentroid'],
                      irafsources['ycentroid'],
                      irafsources['flux']])
    
# =============================================================================
#     flux_array=[]
#     id_array=[]
#     xcenter_array=[]
#     ycenter_array=[]
#     flux_unc_array=[]
# =============================================================================
    
    residual_image=imgdata
    
    

    daogroup = DAOGroup(3 * fwhm)
    psf_model = IntegratedGaussianPRF(sigma=fwhm * gaussian_fwhm_to_sigma)
    psf_model.x_0.fixed = True
    psf_model.y_0.fixed = True
    
    #
    photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=None,
                                    psf_model=psf_model,
                                    fitter=fitter,
                                    fitshape=int(3*fwhm),
                                    aperture_radius=(3*fwhm))
    star_groups=Table(names=['id','group_id','x_0', 'y_0', 'flux_0'],
                data=[irafsources['id'],
                        irafsources['id'],
                            irafsources['xcentroid'],
                              irafsources['ycentroid'],
                              irafsources['flux']])
    
    photometry_result = photometry(image=imgdata-bkg, init_guesses=pos)
# =============================================================================
#     
#     
#     
#     flux_array.append(photometry_result_draft['flux_fit'].value[0])
#     flux_unc_array.append(photometry_result_draft['flux_unc'].value[0])
#     id_array.append(photometry_result_draft['id'].value[0])
#     xcenter_array.append(photometry_result_draft['xcenter'].value[0])
#     ycenter_array.append(photometry_result_draft['ycenter'].value[0])
#     
#     photometry_result=(Table())
#     photometry_result['flux_unc']=flux_unc_array
#     
#     photometry_result['flux_0']=flux_array
#     photometry_result['flux_fit']=flux_array
#     
#     photometry_result['id']=id_array
#     photometry_result['x_0']=xcenter_array
#     photometry_result['y_0']=ycenter_array
#     
#     photometry_result['x_fit']=xcenter_array
#     photometry_result['y_fit']=ycenter_array
#     
#     
#     #Get Residual Image
# =============================================================================
    if produce_residual_data:
        residual_image=photometry.get_residual_image()
    
        fits_reisdual_image=fits.PrimaryHDU(data=residual_image.astype(type=imgdata.dtype),header=hdr)
        fits_reisdual_image.writeto((filepath.split('.fits')[0]+'_2_residual.fits'),overwrite=True)
    
    return photometry_result


def perform_PSF_photometry_sat(sat_x, sat_y, fwhm, imgdata, bkg,
                           fitter=LevMarLSQFitter(), fitshape=5):
    """
    Perform PSF photometry on all sources in a selected image.

    Parameters
    ----------
    irafsources : astropy.table.Table
        Table containing information of all stars detected in the image.
        Has columns:
            id,
            xcentroid,
            ycentroid,
            fwhm,
            sharpness,
            roundness,
            pa,
            npix,
            sky,
            peak,
            flux,
            mag
    fwhm : float
        Mean FWHM of all sources in the image.
    imgdata : numpy.ndarray
        Data from the fits file.
    bkg : float, optional
        Background value of the image in ADU. The default is None.
    bkg_estimator : callable, instance of any
    photutils.background.BackgroundBase subclass, optional
        bkg_estimator should be able to compute either a scalar background or a
        2D background of a given 2D image.
        If None, no background subtraction is performed. The default is None.
    fitter : astropy.modeling.fitting.Fitter instance, optional
        Fitter object used to compute the optimized centroid positions and/or
        flux of the identified sources.
        The default is LevMarLSQFitter().
    fitshape : int or length-2 array-like, optional
        Rectangular shape around the center of a star which will be used to
        collect the data to do the fitting.
        Can be an integer to be the same along both axes. For example, 5 is
        the same as (5, 5),
        which means to fit only at the following relative pixel positions:
            [-2, -1, 0, 1, 2].
        Each element of fitshape must be an odd number. The default is 25.

    Returns
    -------
    photometry_result : astropy.table.Table or None
        Table with the photometry results, i.e., centroids and fluxes
        estimations and the initial estimates used to
        start the fitting process. Uncertainties on the fitted parameters are
        reported as columns called
        <paramname>_unc provided that the fitter object contains a dictionary
        called fit_info with the key param_cov,
        which contains the covariance matrix. If param_cov is not present,
        uncertanties are not reported.

    """
    daogroup = DAOGroup(2 * fwhm)  # The 2 is critical seperation
    psf_model = IntegratedGaussianPRF(sigma=fwhm * gaussian_fwhm_to_sigma)
    psf_model.x_0.fixed = True
    psf_model.y_0.fixed = True
    pos = Table(names=['x_0', 'y_0'],
                data=[sat_x, sat_y])

    photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=None,
                                    psf_model=psf_model,
                                    fitter=fitter,
                                    fitshape=fitshape)
    photometry_result = photometry(image=imgdata - bkg, init_guesses=pos)
    return photometry_result

def perform_aperture_photometry(irafsources,
                                fwhms,
                                imgdata,
                                bkg_2d,
                                filepath,
                                hdr,
                                produce_residual_data=False,
                                calculate_local_error=True,
                                aperture_estimation_mode='mean',
                                **kwargs):
    
    '''
    Performs Aperture Photometry with an option to use either local background 
    or global background values 
    
    Parameters
    ---------
    
    Inputs:
        irafsources: Astropy.Table
                    Contains:
                        - id: 
                            unique object identification number
                        - xcentorid,ycentroid: 
                            estimate of the x and y centroid of the extracted source (star)
                        - fwhm:
                            estimate of the fwhm of the extracted source
                        - sharpness: 
                            sharpness of the source
                        - roudness: 
                            object roundness
                        - flux: 
                            the object instrumental flux
                        - mag: 
                            objects insutrmental magntiude (-2.5*log10(flux))
                        - **kwargs:
                            can include inner and outer radius
            
                    
        imgdata:
            np.array. The img data obtained from opening the filepath
            
        hdr:
            Astropy Header. Header data of the image
        
        bkg 2d: Background2D object
            Describes the backgrounds calculated in previous steps
        
        bkg_std: Global background standard deviation multiplied by the image size. 
        Used to calculate the error in the aperture with the function calc_total_error
        
        produce_residual_data: Boolean
        Defualt is False. Used mainly for debugging purposes to identify which
        stars have not been captured in the image. 
        Residual Data is created by multiplying the Aperture mask by the
        maximum value encountered in the mask and subtracting it from the imgdata. 
        This creates the effect such that all aparture photometry effects can be seen
        when viewing the negative values of the data
        
        calculate_local_error: Boolean
        Aperture Photometry Typically uses the local background to subtract from
        imgdata to get the aperture_sum (flux) amount. Setting boolean to False
        utilizes the global background value calculated in previous steps.

        aperture_estimation_mode: string
        Options are either 'mean' or 'dynamic'
            dynamic mode adjusts the aparture according to the individual fwhm
            mean mode takes the average fwhm and creates an aperture from that
    '''
    

    
    positions=[irafsources['xcentroid'],irafsources['ycentroid']]
    
    #apertures=[CircularAperture[positions,r=r]for r in radii]
    apertures=[]
    #table_result=Table(names=['id','xcenter','ycenter','aperture_sum'],dtype=['int32','float64','float64','float64'])
    
    #table_array=(np.array(aperture_photometry(imgdata-bkg,aperture,error=(bkg_std))))
    #flux_array=[ApertureStats(imgdata-bkg,aperture).std]
    
    flux_unc_array=[]
    flux_array=[]
    id_array=[]
    xcenter_array=[]
    ycenter_array=[]
    masks=np.zeros(np.shape(imgdata),dtype=bool)
    residual_image=imgdata
    
    # TODO: come up with better method for this
    
    # sigma values recommended by Howell
    sigclip= SigmaClip(sigma=3,maxiters=10)
    if 'EGAIN' in hdr:
        gain=hdr['EGAIN']
    elif 'GAINADU' in hdr:
        gain=hdr['GAINADU']

    # TODO: Add  a more dynamic way for initializing the table
    photometry_result=Table(names=['id','xcenter','ycenter','aperture_sum','aperture_sum_err'],
                            dtype=['int','float64','float64','float64','float64'])
    if aperture_estimation_mode == 'manual':
        if 'inner_rad' in kwargs.keys() and 'outer_rad' in kwargs.keys():
            positions_list = [(irafsource['xcentroid'], irafsource[
                'ycentroid']) for irafsource in irafsources]

            apertures= CircularAperture(positions_list,r=kwargs.get('inner_rad'))
            mask1 = apertures.to_mask(method='exact')
            annulus_apertures = CircularAnnulus(positions_list, r_in=kwargs.get('inner_rad'),
                                                r_out=kwargs.get('outer_rad'))
            mask2 = annulus_apertures.to_mask(method='exact')
            sigclip = SigmaClip(sigma=3.0, maxiters=10)
            bkg_stats = ApertureStats(imgdata, annulus_apertures, sigma_clip=sigclip)
            aper_stats = ApertureStats(imgdata, apertures, sigma_clip=None)
            bkg_mean = bkg_stats.median
            aperture_area = apertures.area_overlap(imgdata)
            total_bkg = bkg_mean * aperture_area
            if 'gain' in locals():
                error = calc_total_error(imgdata - bkg_2d.background, bkg_2d.background_rms, gain)
                photometry_result = aperture_photometry(imgdata - bkg_2d.background, apertures, error=error)
            else:
                # Method 1: Calculate standard deviation
                phot_table = aperture_photometry(imgdata-bkg_2d.background, apertures)
                #photometry_result = phot_table
                #photometry_result['aperture_sum'] = phot_table['aperture_sum'] #- total_bkg



                # Calculate error manually, assuming bkg_std is the same
                photometry_result['aperture_sum_err'] = np.sqrt(np.array(apertures.area_overlap(imgdata)) * (
                        bkg_2d.background_rms_median ** 2))


                # Calculate the local aperture sum error manually without effective gain.
                # Aperture Sum Error equation taken from ApertureStats Astropy notebook
                photometry_result['aperture_sum_err']=np.array([np.sqrt(aperture.area_overlap(imgdata)*
                                                                          (np.median(aperture.to_mask().
                                                                                     get_values(bkg_2d.
                                                                                                background_rms))**2))
                                                                                                    for aperture in apertures])




        else:
            print('manual aperture keywords not found, reverting to mean aperture')
            aperture_estimation_mode='mean'

    if aperture_estimation_mode == 'mean':

        fwhm=np.mean(fwhms)
        positions_list=[(irafsource['xcentroid'],irafsource[
            'ycentroid']) for irafsource in irafsources]
        apertures=CircularAperture(positions_list,
                                   r=3*fwhm)
        mask1=apertures.to_mask(method='exact')


        annulus_apertures=CircularAnnulus(positions_list,r_in=3*fwhm,
                                         r_out=4*fwhm)
        mask2=annulus_apertures.to_mask(method='exact')
        sigclip=SigmaClip(sigma=3.0,maxiters=10)

        bkg_stats=ApertureStats(imgdata,annulus_apertures,sigma_clip=sigclip)
        aper_stats=ApertureStats(imgdata,apertures,sigma_clip=None)
        bkg_mean=bkg_stats.median
        aperture_area=apertures.area_overlap(imgdata)
        total_bkg=bkg_mean*aperture_area


        ### Method One: Take the Standard deviation of the aperture itself
        # phot_table = aperture_photometry(imgdata, apertures)
        # photometry_result = phot_table
        # photometry_result['aperture_sum'] = phot_table[
        #                                        'aperture_sum'] - total_bkg
        # photometry_result['aperture_sum_err'] = aper_stats.std

        #Work around way to calculate the total error using localized background standard deviation

        ### Method 2: Calculate Total Error

        if 'gain' in locals():
            # FIXME : Add Local Error Code for this condition
            error = calc_total_error(imgdata - bkg_2d.background, bkg_2d.background_rms, gain)
            photometry_result = aperture_photometry(imgdata - bkg_2d.background, apertures, error=error)
        else:
            # Method 1: Calculate standard deviation
            phot_table = aperture_photometry(imgdata-bkg_2d.background, apertures)
            photometry_result = phot_table
            photometry_result['aperture_sum'] = phot_table['aperture_sum'] #- total_bkg

            #Calculate error manually, assuming bkg_std is the same throughout the image
            #TODO : Double check the calculations for this method
            photometry_result['aperture_sum_err'] = np.sqrt(np.array(apertures.area_overlap(imgdata))*(
                    bkg_2d.background_rms_median**2))
            # Aperture Sum Error equation taken from ApertureStats Astropy notebook

        if produce_residual_data:
            masks1= [mask.to_image(np.shape(imgdata)) for mask in mask1]
            masks2 =[mask.to_image(np.shape(imgdata)) for mask in mask2]
            masks = (np.sum(masks1,axis=0) + np.sum(masks2,axis=0))

            # Creates a rough reisudal imae to show where the aperture masks are
            residual_image = (imgdata - (
            (np.max(imgdata) * masks))).astype(
                'float32')



    if aperture_estimation_mode == 'dynamic':
        if 'gain' in locals():
            error = calc_total_error(imgdata - bkg_2d.background, bkg_2d.background_rms, gain)
        for i in range(0,np.shape(positions)[1]):
            radii=3*fwhms[i]

            aperture=CircularAperture((positions[0][i],positions[1][i]),r=radii)
            apertures.append(aperture)


            mask1=aperture.to_mask(method='exact')

            combined_aperture_mask=mask1.to_image(np.shape(imgdata))


            # Test Annulus inner and Outer dimensions  See Shaw pg. 55
            "Calculate Local Error"
            annulus_aperture=CircularAnnulus((positions[0][i],positions[1][i]),r_in=3*radii,r_out=4*radii)

            #bkg_stats=ApertureStats(imgdata,annulus_aperture,sigma_clip=sigclip)
            #bkg = bkg_stats.median
            #aper_stats = ApertureStats(imgdata-bkg, aperture)

            if produce_residual_data:
                # Masking
                mask2=annulus_aperture.to_mask(method='exact')
                combined_aperture_mask=combined_aperture_mask+mask2.to_image(np.shape(imgdata))






            ### Method 2: Calculate Total Error
            # bkg_error = (np.ones(np.shape(imgdata))) * bkg_stats.std
            if 'gain' in locals():

                photometry_result_draft=aperture_photometry(imgdata-bkg_2d.background,aperture,error=error)

            elif 'gain' not in locals(): # Poission gain not taken into effect when gain is not in FITS Headers
                photometry_result_draft = aperture_photometry(imgdata - bkg_2d.background, aperture)

                # Calculate the local aperture sum error manually without effective gain.
                # Aperture Sum Error equation taken from ApertureStats Astropy notebook
                photometry_result['aperture_sum_err'] = np.sqrt(aperture.area_overlap(imgdata) *
                                                                          (np.median(aperture.to_mask().
                                                                                     get_values(bkg_2d.
                                                                                                background_rms)) ** 2))
            photometry_result_draft[0]['id']=i+1
            photometry_result.add_row(photometry_result_draft[0])


            if produce_residual_data:

                "Masks Used for Debugging"

                masks=(masks+combined_aperture_mask)

            # Creates a rough residual imae to show where the aperture masks are
                residual_image=(imgdata-((np.max(imgdata)*
                        combined_aperture_mask))).astype('float32')

    if produce_residual_data:
        "Mask and Residual Image Output"
        hdr['IMAGETYP']='MASK'
        mask_image=fits.PrimaryHDU(data=masks.astype(dtype=imgdata.dtype),header=hdr)
        mask_image.writeto((filepath.split('.fits')[0]+'_mask.fits'),overwrite=True)
        hdr['IMAGETYP']='RESIDUAL'
        residual_image=fits.PrimaryHDU(data=residual_image,header=hdr)
        residual_image.writeto((filepath.split('.fits')[0]+'_aperture_residual.fits'),overwrite=True)

    return photometry_result
