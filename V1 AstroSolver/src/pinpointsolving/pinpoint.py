

import os

import win32com
import win32com.client

import AstroFunctions as astro


def pinpoint_solve(inbox, catloc, max_mag=12, sigma=5.0, catexp=0.4, match_residual=1.5,
                   max_solve_time=300, space_based_bool=False, use_sextractor=False,
                   all_sky_solve=False):

    file_suffix = (".fits", ".fit", ".fts",".FIT",".FITS")
    failedSolves=0
    for dirpath, dirnames, filenames in os.walk(inbox):
        for filename in filenames:
            if (filename.endswith(file_suffix)):
                filepath = os.path.join(dirpath, filename)
                # Creates an instance of a plate object
                f = win32com.client.Dispatch("Pinpoint.plate")
                print("Processing Image: " + filepath)

                "Import Data from FITS Image"
                __, __, __, __, __,\
                    __, __, header, __, __, __ = \
                    astro.fits_header_import(filepath, space_based_bool)

                """Pinpoint Solve"""
                try:
                    # Attaches FITS image to the created Plate
                    f.AttachFITS(filepath)
                    f.Declination = f.targetDeclination
                    f.RightAscension = f.targetRightAscension
                    try:
                        x_arcsecperpixel, y_arcsecperpixel =\
                            astro.calc_ArcsecPerPixel(header)
                        # f.ArcsecperPixelHoriz = x_arcsecperpixel
                        # f.ArcsecperPixelVert = y_arcsecperpixel
                    except:

                        # yBin = 4.33562092816E-004*3600;
                        # xBin =  4.33131246330E-004*3600;
                        f.ArcsecperPixelHoriz = 1.38
                        f.ArcsecperPixelVert = 1.38


                    "Pinpoint Solve Inputs"
                    f.Catalog = 11
                    f.CatalogPath = catloc
                    f.UseSExtractor = use_sextractor
                    f.CatalogMaximumMagnitude = max_mag
                    f.CatalogExpansion = catexp
                    f.MaxSolveTime = max_solve_time
                    f.MaxMatchResidual = match_residual
                    f.SigmaAboveMean = sigma
                    f.RemoveHotPixels()
                    f.Solve()
                    f.UpdateFITS()

                    print("Solved - Header Updated")
                except:
                    print("Could Not Solve")
                    failedSolves+=1


    return None
def pinpoint_init():
    f = win32com.client.Dispatch("Pinpoint.plate")
    return f