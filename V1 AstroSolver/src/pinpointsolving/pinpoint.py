

import os

import win32com
import win32com.client

import AstroFunctions as astro


def pinpoint_solve(inbox, catloc, max_mag=13, sigma=2.5, catexp=0.8, match_residual=0,
                   max_solve_time=300, space_based_bool=False, use_sextractor=False,
                   all_sky_solve=False):

    file_suffix = (".fits", ".fit", ".fts")

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

                    x_arcsecperpixel, y_arcsecperpixel =\
                        astro.calc_ArcsecPerPixel(header)
                    # yBin = 4.33562092816E-004*3600;
                    # xBin =  4.33131246330E-004*3600;
                    # f.ArcsecperPixelHoriz  = 4.556
                    # f.ArcsecperPixelVert =  4.556
                    f.ArcsecperPixelHoriz = x_arcsecperpixel
                    f.ArcsecperPixelVert = y_arcsecperpixel

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

                    print("Pinpoint Solved")
                except:
                    print("Could Not Solve")
                    continue
    return f
def pinpoint_init():
    f = win32com.client.Dispatch("Pinpoint.plate")
    return f