#Subtract background first
#Compute the Radon Transform using FFT
    #Compute the length and orientation of the streaks
#Filter star streaks with an iterative matched filter
#Non-streak objects (PSFs) detected by subtracting the filtered streak only
    #  image from the background-free image


import os
import numpy
import numpy as np
import scipy.signal
import skimage
from skimage.data import shepp_logan_phantom
from skimage.transform import radon, rescale
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft
import scipy.fft
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
streak=r'/Users/home/Downloads/April21/TRM'
file_suffix = (".fits", ".fit", ".FIT")

def ImageSize(data):
    x,y=numpy.size(data)
    return x,y
def killBrightSpots(data, sigma=3.0):
    '''
    Remove all pixels with a brightness greater than sigma*sigma*mean(data)
    Parameters
    ----------
    data
    sigma

    Returns
    -------

    '''

    x,y=ImageSize(data)
    spots=numpy.zeros(x,y)
    imgout=data
    outwindowindex=[-7,-4,  -1,  1,  4, 7]
    k3 = numpy.ones(3, 3)
    k3 = k3/9
    noise=sigma*10
    means=scipy.signal.convolve2d(data,k3)
    for i in range(9,x-8):
        for j in range(9,y-8):
            if means[i,j]>noise:
                out=means[i+outwindowindex,j+outwindowindex]
                out[2:5,2:5]=0
                maxmeanout=max(numpy.mean(out))
                if (5*maxmeanout)<means[i,j]:
                    spots[i-3:i+3,j-3:j+3]=data[i-3:i+3,j-3:j+3]
                    imgout[i-3:i+3,j-3:j+3]=0
    return imgout, spots

def background_estimation(data):
    return

def fix_dropped_pixels(data,sigma):
    dataOutput=data
    droppedrows=[]
    x,y=ImageSize(data)
    maskx, masky = numpy.where(data==0)
    im_as_row=numpy.reshape(data,1,x*y)
    a,b=numpy.histogram((im_as_row,numpy.unique(im_as_row)))
    mode=b(a==max(a))
    if mode==0:
        zeroIndex=numpy.where(a==max(a))
        a[zeroIndex]=0
        mode=b(a==max(a))
    mode=max(mode)
    rowSum=sum(numpy.array(data).H)
    for i in range(1,x):
        if not rowSum[i]:
            droppedrows=[droppedrows,i]
    for i in range(1,ImageSize(maskx)):
        dataOutput[maskx[i],masky[i]]=mode
    return dataOutput, droppedrows

def smoothing(A,steps,kernalsize):
    B=A
    Asmoothed=B
    N,y=max(ImageSize(A))
    for i in range(1,steps):
        if kernalsize==2:
            Asmoothed[2:N-1]=B[1:N-2]/4+B[2:N-1]/2+B[3:N]/4
            Asmoothed[1]=B[1]
            Asmoothed[N]=B[N]
        elif kernalsize==3:
            Asmoothed[3:N-2] = B[1:N-4]/9 + B[2:N-3]/4.5 + B[3:N-2]/2
            Asmoothed[1] = B[1]
            Asmoothed[2] = B[1]/4 + B[2]/2 + B[3]/4
            Asmoothed[N-1] = B[N-1]/2 + B[N-1]+ B[N]/4
            Asmoothed[N] = B[N]
        else:
            print("Kernal Size above 3")
        B=Asmoothed
    return B

def find_sinxx2_maxima(sinxx,origin_index,meanfreq,halfband, N):
    indl=numpy.floor(origin_index + (N*meanfreq)-numpy.floor(halfband))
    indr=numpy.ceil(origin_index+(N*meanfreq)+numpy.floor(halfband))
    localmax=max(sinxx[indl:indr])
    index=indl-1 + numpy.where(sinxx[indl:indr]==localmax)
    if index>indl and index<indr:
        freq=index-origin_index
    else:
        index = origin_index
        freq=0
    return index, freq

def find_sinxx2_minima(sinxx,origin_index,meanfreq,halfband, N):
    indl = numpy.floor(origin_index + (N*meanfreq) - numpy.floor(halfband))
    indr = numpy.ceil (origin_index + (N*meanfreq) + numpy.floor(halfband))
    localmin = min(sinxx[indl:indr])
    index = indl-1 + numpy.where(sinxx[indl:indr]==localmin)
    if index>indl and index<indr:
        freq=index-origin_index
    else:
        index = origin_index
        freq=0
    return index, freq

def quickradonsearchiter(fftabs,R,rlindex,anglemax,rmax,angle_increment,rprofile):
    anglemid=anglemax
    maxmid=rmax

    angleleft=numpy.mod(anglemid - angle_increment,180)
    angleright = numpy.mod(anglemid + angle_increment, 180)

    radonleft=radon(fftabs,angleleft)
    radonright=radon(fftabs,angleright)

    maxleft=radonleft[rlindex]
    maxright=radonright[rlindex]

    angles=[angleleft,anglemid,angleright]
    maxes=[maxleft,maxmid,maxright]

    rmax=max(maxes)
    max_index=numpy.where(maxes=rmax)
    anglemax=angles[max_index]
    profiles=[radonleft,rprofile,radonright]
    rprofile=profiles[:max_index]
    return rprofile

def quickradonsearch(fftabs):
    R=radon(fftabs,theta=[0,20,40,60,80,100,120,140,160,180])
    rx,ry=ImageSize(R)
    rlindex=numpy.floor((rx+1)/2)
    rmax=max(R[rlindex][:])
    angleIndex=numpy.where(R[rlindex][:]==rmax)
    anglemax=(angleIndex-1)*20
    rprofile=R[:][angleIndex]
    angle_increment=10
    while angle_increment>0.1:
        R=quickradonsearchiter()


def find_sinxx2_extreme(rprofile):
    sinxx=smoothing(rprofile,2,3)
    sinxx2=sinxx
    nmax7_frequency = 0
    nmax5_frequency = 0
    nmax3_frequency = 0
    pmax3_frequency = 0
    pmax5_frequency = 0
    pmax7_frequency = 0
    nmin6_frequency = 0
    nmin4_frequency = 0
    nmin2_frequency = 0
    pmin2_frequency = 0
    pmin4_frequency = 0
    pmin6_frequency = 0
    originVal=max(sinxx)
    originIndex=numpy.where(sinxx==originVal)
    if ImageSize(originIndex,1)>1: #TODO add multi dimension option to ImageSize
        #TODO add line 281 error catch
    clip_profile= sinxx>=originVal/2
    hits=numpy.where(clip_profile==1)
    bumps=numpy.zeros(ImageSize(hits))
    for i in range(2,ImageSize(hits)):
        bumps[i]=hits[i]-hits[i-1]
    cuts=numpy.where(bumps>1)
    if cuts and (ImageSize(hits)>ImageSize(clip_profile)/2):
        ncuts=ImageSize(cuts)
        first_cut=cuts[1]
        last_cut=cuts[-1]
        for i in range(1,first_cut-1):
            clip_profile[hits[i]]=0
        for i in range(last_cut,ImageSize(hits)):
            clip_profile[hits[i]]=0

        originVal=max(sinxx*clip_profile)
        originIndex=numpy.where(sinxx==originVal)
        if ImageSize(originIndex)>1:
            originIndex=round(numpy.mean(originIndex))
    hlhw_left=numpy.where(clip_profile, 1, 'first')
    hlhw_right=numpy.where(clip_profile, 1, 'last')
    hlhw=hlhw_right-hlhw_left

    lsxx = len(sinxx)
    sinxxh = sinxx* hanning(lsxx)
    Y = fft(sinxxh)
    Mag = abs(Y[1:numpy.floor(lsxx / 2)])**2
    b = numpy.where(diff(Mag) > 0, 1)
    if hlhw > (6 * (8 / 3) * b):
        hlhw = (6 * (8 / 3) * b)



def radon_analysis(fftabs):
    x,y=ImageSize(fftabs)
    
def make_kernal_line(angle,length,option):

    angle=numpy.arctan(numpy.tan(angle))
    dx=numpy.ceil(length*max(np.abs(np.sin(angle)),np.cos(angle)))
    dy=dx
    linekernal = np.zeros(dx,dy)
    cx=(dx+1)/2
    cy=cx

    if ((angle > (np.pi/4)) or (angle < (-np.pi/4))):

        vertIncrement = np.cos(angle)/np.sin(angle)
        for j in range(dy):
            x = np.round(cx+(j-cy)*vertIncrement)
            linekernal[x][j] = 1
    else:
        horizIncrement = np.sin(angle) / np.cos(angle)
        for i in range(dx):
            x = np.round(cy + (j - cx) * horizIncrement)
            linekernal[i][y] = 1

    linekernal = linekernal/np.sum(np.sum(linekernal))
    return linekernal


def make_matched_filter(kernal, xsize, ysize):
    row,col=ImageSize(kernal)
    midrow = np.floor(row/2)
    midcol = np.floor(col/2)
    if xsize>row and ysize>col:
        filter = np.zeros()
        filter[1:row-midrow][1:col-midcol] = kernal[midrow+1:row][midcol+1:col]
        filter[1:row - midrow][ysize-midcol+1:ysize] = kernal[midrow + 1:row][1:col]
        filter[xsize - midrow+1:xsize][ysize-midcol+1:ysize] = kernal[1:midrow][1:midcol]
        filter[xsize - midrow+1:xsize][1:col - midcol] = kernal[1:midrow][midcol + 1:col]
    ft=fft.fft2(filter)
    return ft

def match_streak(image,streakAngle,streakLength,niter):
    kernal=make_kernal_line(streakAngle, streakLength)
    filter=make_matched_filter(kernal, ImageSize(image))
    clippedImage=image
    for i in range(niter):
        ftimage=scipy.fft2(clippedImage)
        convolve=np.abs(scipy.ifft(ftimage*filter))
        clippedImage=min(image,2*convolve)
    return clippedImage, convolve

def match_streak2(image,streakAngle,streakLength,niter):
    kernal=make_kernal_line(streakAngle, streakLength)
    filter=make_matched_filter(kernal, ImageSize(image))
    clippedImage=image
    for i in range(niter):
        convolved_image=convolve_fft(clippedImage,filter)
        clippedImage=min(image, convolved_image)
    return clippedImage, convolved_image



def measureStreaks(image, streakAngle, streakLength):
    newImage, convolvedImage = match_streak(streakAngle,streakLength)


for dirpath, dirnames, filenames in os.walk(streak):
    print(f'{len(filenames)} Images Detected')
    for filename in filenames:
        if (filename.endswith(file_suffix)):
            # sum1+=1
            filepath = os.path.join(dirpath, filename)
            STARS = open(filepath + '.stars', "w")
            AstrometryNetFile = open(filepath + '.txt', 'w')

            imagehdularray = fits.open(filepath)
            fitsdata = imagehdularray[0].data
            sigma_clip = SigmaClip(sigma=3)
            try:
                bkg = Background2D(fitsdata, (30, 30), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
            except:
                bkg = SExtractorBackground(fitsdata)
            nobkgfits= fitsdata - bkg.background

            #Fourier Transform the Image


            #Radon Transform
            x,y=numpy.size(fftabs)
            data,spots=killBrightSpots(nobkgfits)
            fft=numpy.ifftshift(nobkgfits)
            fft = numpy.fft2(fft)
            fftshift = numpy.fftshift(fft)
            fftabs = numpy.abs(fftshift)

