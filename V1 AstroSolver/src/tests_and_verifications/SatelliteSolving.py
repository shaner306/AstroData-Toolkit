import os
import numpy
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch
from astropy.visualization import simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.background import Background2D
from photutils.background import SExtractorBackground
from photutils.segmentation import SourceCatalog
from photutils.segmentation import deblend_sources
from photutils.segmentation import detect_sources
from photutils.segmentation import detect_threshold
from photutils.segmentation import make_source_mask

import StreakDetectionTest
import pinpointsolving.pinpoint

streak =r'/Users/home/Sync/'
ref_stars_file = r'D:\Astro2\Reference Star Files\Reference_stars_Apr29.txt'
catloc1 = "D:\\squid\\USNOA20-All"
catloc = 'D:\\squid\\UCAC4'


StreakDetectionTest.StreakDetection(streak, TRM=True)

pinpointsolving.pinpoint_solve(streak, catloc, max_mag=12, sigma=5.0, catexp=0.4, match_residual=2.5,
                               max_solve_time=300, space_based_bool=False, use_sextractor=True,
                               all_sky_solve=False)



