

import os

import win32com
import win32com.client
import AstroFunctions as astro


def pinpoint_solve(inbox, catloc, max_mag, sigma, catexp, match_residual,
                   max_solve_time, cat, space_based_bool, use_sextractor,
                   all_sky_solve,pinpoint=True):

    file_suffix = (".fits", ".fit", ".fts")

    for dirpath, dirnames, filenames in os.walk(inbox):
        for filename in filenames:
            if (filename.endswith(file_suffix)):

                filepath = os.path.join(dirpath, filename)
                # Creates an instance of a plate object
                f = win32com.client.Dispatch("Pinpoint.plate")

                print("Processing Image: " + filepath)

                "Import Data from FITS Image"
                header = 0
                imagehdularray, date, exposure_Time, imagesizeX, imagesizeY,\
                    fitsdata, filt, header, XPIXSZ, YPIXSZ, wcs = \
                    astro.fits_header_import(filepath, space_based_bool)

                """Pinpoint Solve"""
                # FIXME: Why is pinpont variable here when always TRUE?
                if pinpoint:
                    try:
                        # Attaches FITS image to the created Plate
                        f.AttachFITS(filepath)
# We set the declination of the image to the target object's declination
                        f.Declination = f.targetDeclination
                        f.RightAscension = f.targetRightAscension

                        x_arcsecperpixel, y_arcsecperpixel =\
                            astro.calc_ArcsecPerPixel(header)
                        # yBin = 4.33562092816E-004*3600;
                        # xBin =  4.33131246330E-004*3600;
                        # f.ArcsecperPixelHoriz  = 4.556
                        # f.ArcsecperPixelVert =  4.556
                        # if space_based_bool==1:
                        #     yBin = 4.33562092816E-004*3600*2
# %Image Specific Pixel Size in arcsec /
# Obtained from FITS Header and converted from deg
                        #     xBin =\
# 4.33131246330E-004*3600*2#%Image Specific Pixel Size in arcsec /
# Obtained from FITS Header and converted from deg
                        # else:
                        #     yBin = 4.33562092816E-004*3600
                        # %Image Specific Pixel Size in arcsec /
                        # Obtained from FITS Header and converted from deg
                        #     xBin =  4.33131246330E-004*3600
                        # %Image Specific Pixel Size in arcsec /
                        # Obtained from FITS Header and converted from deg
                        # f.ArcsecperPixelHoriz  = xBin
                        # %CCD Pixel scale on CD
                        # f.ArcsecperPixelVert = yBin
                        f.ArcsecperPixelHoriz = x_arcsecperpixel
                        # %CCD Pixel scale on CD
                        f.ArcsecperPixelVert = y_arcsecperpixel

                        "Pinpoint Solve Inputs"
                        # TODO Add Inputs for pinpoint solving to GUI
                        f.Catalog = 11
                        f.CatalogPath = catloc
                        f.UseSExtractor = use_sextractor
                        f.CatalogMaximumMagnitude = max_mag
                        # print(f.CatalogMaximumMagnitude)
                        f.CatalogExpansion = catexp
                        f.MaxSolveTime = max_solve_time
                        f.MaxMatchResidual = match_residual
                        f.SigmaAboveMean = sigma
                        f.RemoveHotPixels()

                        "Pinpoint Solving"

                        # FIXME f.Solve intronsicly calls Find Catalog Stars and Image Stars
                        # f.FindCatalogStars()
                        # print(f.CatalogStars.Count)
                        # f.FindImageStars()
                        # print(f.ImageStars.Count)
                        f.Solve()
                        f.UpdateFITS()
                        # print( f.MatchedStars.count)
                        # f.FindImageStars()
                        # print(f.ImageStars)

                        f.DetachFITS()
                        f = None
                        print("Pinpoint Solved")
                    except:
                        print("Could Not Solve")
                        continue