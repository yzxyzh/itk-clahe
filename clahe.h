/**
  \file clahe.h
  \date 8 april 2011
  \author Francis Girard
  
  Copyright (C) 2010 Francis Girard
  
  GNU LESSER GENERAL PUBLIC LICENSE
  Version 3, 29 June 2007
  
  Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
  
  Everyone is permitted to copy and distribute verbatim copies of this license
  document, but changing it is not allowed.
  
  This version of the GNU Lesser General Public License incorporates the terms
  and conditions of version 3 of the GNU General Public License, supplemented by
  the additional permissions listed below.
  
  The precise terms and conditions for copying, distribution and modification
  may be found here:
  
  http://www.gnu.org/licenses/lgpl.html
  
  
  
  Contrast Limited Adaptive Histogram Equalization
  
  Adapted from original Karel Zuiderveld's code for ITK with minor
  optimizations here and there.
  
  Original C file is without licence and contains these heading comments:
  
  ANSI C code from the article
  "Contrast Limited Adaptive Histogram Equalization"
  by Karel Zuiderveld, karel@cv.ruu.nl
  in "Graphics Gems IV", Academic Press, 1994
  
  These functions implement Contrast Limited Adaptive Histogram Equalization.
  The main routine (CLAHE) expects an input image that is stored contiguously in
  memory;  the CLAHE output image overwrites the original input image and has the
  same minimum and maximum values (which must be provided by the user).
  This implementation assumes that the X- and Y image resolutions are an integer
  multiple of the X- and Y sizes of the contextual regions. A check on various other
  error conditions is performed.
  
  #define the symbol BYTE_IMAGE to make this implementation suitable for
  8-bit images. The maximum number of contextual regions can be redefined
  by changing uiMAX_REG_X and/or uiMAX_REG_Y; the use of more than 256
  contextual regions is not recommended.
  
  The code is ANSI-C and is also C++ compliant.
  
  Author: Karel Zuiderveld, Computer Vision Research Group,
          Utrecht, The Netherlands (karel@cv.ruu.nl)

*/

#ifndef __clahe_h__
#define __clahe_h__

#include <math.h>
#include "itkImage.h"
#include "itkImportImageFilter.h"


/**
  Implementation of the CLAHE algorithm as described by Karel Zuiderveld.
*/
template<typename T_Pixel>
class ClaheImpl
{
  public:
    
    inline ClaheImpl();
    
    inline int execute
    (
      T_Pixel* pImage,
      unsigned int uiXRes,
      unsigned int uiYRes,
      T_Pixel Min,
      T_Pixel Max,
      unsigned int uiNrX,
      unsigned int uiNrY,
      unsigned int uiNrBins,
      float fCliplimit
    );
    
    inline int wrapedExecute
    (
      T_Pixel* pImage,
      unsigned int uiXRes,
      unsigned int uiYRes,
      T_Pixel Min,
      T_Pixel Max,
      unsigned int uiNrX,
      unsigned int uiNrY,
      unsigned int uiNrBins,
      float fCliplimit
    );
    
  private:
    
    inline void _clipHistogram
    (
      unsigned long* pulHistogram,
      unsigned int uiNrGreylevels,
      unsigned long ulClipLimit
    );
    
    inline void _makeHistogram
    (
      T_Pixel* pImage,
      unsigned int uiXRes,
      unsigned int uiSizeX,
      unsigned int uiSizeY,
      unsigned long* pulHistogram,
      unsigned int uiNrGreylevels,
      unsigned int* pLookupTable
    );
    
    inline void _mapHistogram
    (
      unsigned long* pulHistogram,
      T_Pixel Min,
      T_Pixel Max,
      unsigned int uiNrGreylevels,
      unsigned long ulNrOfPixels
    );
    
    inline void _makeLut
    (
      unsigned int* pLUT,
      T_Pixel nMin,
      T_Pixel nMax,
      unsigned int uiNrBins
    );
    
    inline void _interpolate
    (
      T_Pixel* pImage,
      unsigned int uiXRes,
      unsigned long* pulMapLU,
      unsigned long* pulMapRU,
      unsigned long* pulMapLB,
      unsigned long* pulMapRB,
      unsigned int uiXSize,
      unsigned int uiYSize,
      unsigned int* pLUT
    );
    
    inline unsigned int _pow
    (
      unsigned int nBase,
      unsigned int nExp
    );
    
  private:
    
    /* max. # contextual regions in x-direction */
    static const unsigned int uiMAX_REG_X = 16;
    
    /* max. # contextual regions in y-direction */
    static const unsigned int uiMAX_REG_Y = 16;
};



/**
  Applies Clahe to an ITK image to produce another image.
*/
template<typename T_ItkImage>
class ClaheITK
{
  public:
    
    typedef typename T_ItkImage::PixelType T_Pixel;
    typedef typename T_ItkImage::Pointer T_ItkImagePointer;
    typedef typename T_ItkImage::RegionType T_RegionType;
    typedef typename T_ItkImage::SizeType T_SizeType;
    
    inline ClaheITK();
    inline ~ClaheITK();
    
    inline void SetInput(T_ItkImagePointer pInput);
    inline typename T_ItkImage::Pointer GetOutput();
    inline void Update();
    
    inline void setGrayLevelMin(T_Pixel nMin);
    inline void setGrayLevelMax(T_Pixel nMax);
    inline void setNbRegionsX(unsigned int uiNrX);
    inline void setNbRegionsY(unsigned int uiNrY);
    inline void setNbBins(unsigned int uiNrBins);
    inline void setCliplimit(float fCliplimit);
    inline void setAutoAdaptToNbRegions(bool bAutoAdaptToNbRegions);
    
    inline T_Pixel getGrayLevelMin() const;
    inline T_Pixel getGrayLevelMax() const;
    inline unsigned int getNbRegionsX() const;
    inline unsigned int getNbRegionsY() const;
    inline unsigned int getNbBins() const;
    inline float getCliplimit() const;
    inline bool getAutoAdaptToNbRegions() const;
    
    inline void execute(T_ItkImagePointer pImage);
  
  private:
  
    inline T_Pixel* _fromItkToFlatArray
    (
      T_ItkImagePointer pImage,
      const T_RegionType& roRegion,
      unsigned int nNbBytes
    );
    
    inline void _fromFlatArrayToItk
    (
      T_Pixel* pFlatImage,
      const T_SizeType& roImgSize,
      unsigned int nNbPixels
    );
  
  private:
  
    typedef itk::ImportImageFilter<T_Pixel, 2> T_ImportFilter;
  
  private:
    
    T_Pixel _nMin;
    T_Pixel _nMax;
    unsigned int _uiNrX;
    unsigned int _uiNrY;
    unsigned int _uiNrBins;
    float _fCliplimit;
    bool _bAutoAdaptToNbRegions;
    T_ItkImagePointer _pInput;
    
    T_Pixel* _pFlatImage;
    ClaheImpl<T_Pixel> _oClaheImpl;
    typename T_ImportFilter::Pointer _pImportFilter;
};



/**
  Useless ctor!
*/
template<typename T_Pixel>
inline ClaheImpl<T_Pixel>::ClaheImpl()
{
}


/**
*/
template<typename T_Pixel>
inline int ClaheImpl<T_Pixel>::wrapedExecute
(
  T_Pixel* pImage,
  unsigned int uiXRes,
  unsigned int uiYRes,
  T_Pixel Min,
  T_Pixel Max,
  unsigned int uiNrX,
  unsigned int uiNrY,
  unsigned int uiNrBins,
  float fCliplimit
)
{
  //
  // -1- Compute lines and columns to add if necessary
  //
  // If the desired number of regions according to X or Y is not a multiple
  // of complete image width or height then reallocate the image with extra
  // lines and/or columns with intensity 0.
  //
  bool bRealloc = false;
  
  unsigned int nWidthMissing = uiNrX - uiXRes % uiNrX;
  unsigned int nLeftToInsert = 0;
  unsigned int nRightToInsert = 0;
  
  if (nWidthMissing != 0)
  {
    bRealloc = true;
    nLeftToInsert = nWidthMissing / 2;
    nRightToInsert = nWidthMissing - nLeftToInsert;
  }
  
  unsigned int nHeightMissing = uiNrY - uiYRes % uiNrY;
  unsigned int nUpToInsert = 0;
  unsigned int nBottomToInsert = 0;
  
  if (nHeightMissing != 0)
  {
    bRealloc = true;
    nUpToInsert = nHeightMissing / 2;
    nBottomToInsert = nHeightMissing - nUpToInsert;
  }
  
  unsigned int nNewWidth = uiXRes + nWidthMissing;
  unsigned int nNewHeight = uiYRes + nHeightMissing;
  T_Pixel* pImageNew = pImage;
  
  //
  // -2- Reallocate image if necessary
  //
  if (bRealloc)
  {
    pImageNew = (T_Pixel*)  malloc(sizeof(T_Pixel) * nNewWidth * nNewHeight);
    T_Pixel* pImageNewIt = pImageNew;
    
    unsigned int nIdx = 0;
    unsigned int nSubIdx = 0;
    
    //
    // Add the upper lines if necessary
    //
    for ( ; nIdx < nUpToInsert; nIdx++)
    {
      nSubIdx = 0;
      for ( ; nSubIdx < nNewWidth;  nSubIdx++, pImageNewIt++)
      {
        *pImageNewIt = 0;
      }
    }
    
    //
    // Copy the lines from the original image, adding extra columns if necessary
    //
    T_Pixel* pImageIt = pImage;
    nIdx = 0;
    nSubIdx = 0;
    for ( ; nIdx < uiYRes; nIdx++)
    {
      // Add some column at the beginning if necessary
      nSubIdx = 0;
      for ( ; nSubIdx < nLeftToInsert;  nSubIdx++, pImageNewIt++)
      {
        *pImageNewIt = 0;
      }
      
      // Copy from original
      nSubIdx = 0;
      for ( ; nSubIdx < uiXRes;  nSubIdx++, pImageNewIt++, pImageIt++)
      {
        *pImageNewIt = *pImageIt;
      }
      
      // Add some column at the end if necessary
      nSubIdx = 0;
      for ( ; nSubIdx < nRightToInsert;  nSubIdx++, pImageNewIt++)
      {
        *pImageNewIt = 0;
      }
    }
    
    //
    // Add the lower lines if necessary
    //
    nIdx = 0;
    nSubIdx = 0;
    for ( ; nIdx < nBottomToInsert; nIdx++)
    {
      nSubIdx = 0;
      for ( ; nSubIdx < nNewWidth;  nSubIdx++, pImageNewIt++)
      {
        *pImageNewIt = 0;
      }
    }
    
  }
  
  
  //
  // -3- Apply CLAHE algorithm
  //
  int nReturnCode = execute
                    (
                      pImageNew,
                      nNewWidth,
                      nNewHeight,
                      Min,
                      Max,
                      uiNrX,
                      uiNrY,
                      uiNrBins,
                      fCliplimit
                    );
  
  
  //
  // 4- Copy back to the original image if necessary.
  //
  if (nReturnCode == 0 && bRealloc)
  {
    T_Pixel* pImageNewIt = pImageNew;
    T_Pixel* pImageIt = pImage;
    
    //
    // Skip the upper lines that we did insert.
    //
    pImageNewIt += (nNewWidth * nUpToInsert);
    
    //
    // Copy the original lines, skipping inserted pixels at left and at right
    //
    unsigned int nIdx = 0;
    unsigned int nSubIdx = 0;
    for ( ; nIdx < uiYRes; nIdx++)
    {
      // Skip inserted column at the left if necessary
      pImageNewIt += nLeftToInsert;
      
      // Copy to original
      nSubIdx = 0;
      for ( ; nSubIdx < uiXRes;  nSubIdx++, pImageNewIt++, pImageIt++)
      {
        *pImageIt = *pImageNewIt;
      }
      
      // Skip inserted column at the right if necessary
      pImageNewIt += nRightToInsert;
    }
    
    //
    // Skip the lower lines that we did insert.
    // (Useless and therefore commented out).
    //
    // pImageNewIt += (nNewWidth * nBottomToInsert);
  }
  
  return nReturnCode;
}



/**
  These functions implement Contrast Limited Adaptive Histogram Equalization.
  The main routine (CLAHE) expects an input image that is stored contiguously in
  memory;  the CLAHE output image overwrites the original input image and has the
  same minimum and maximum values (which must be provided by the user).
  This implementation assumes that the X- and Y image resolutions are an integer
  multiple of the X- and Y sizes of the contextual regions. A check on various other
  error conditions is performed.
  
  The number of "effective" greylevels in the output image is set by uiNrBins;
  selecting a small value (eg. 128) speeds up processing and still produce an
  output image of good quality. The output image will have the same minimum and
  maximum value as the input image. A clip limit smaller than 1 results in
  standard (non-contrast limited) AHE.
  
  \param pImage The one channel 2D image as one flat array. The rows are layed
    out one after the other.
  \param uiXRes Complete width of the image.
  \param uiYRes Complete height of the image.
  \param Min Minimum greyvalue of input image (also becomes minimum of output image)
  \param Max Maximum greyvalue of input image (also becomes maximum of output image)
  \param uiNrX Number of contextual regions in the X direction (min 2, max uiMAX_REG_X)
  \param uiNrY Number of contextual regions in the Y direction (min 2, max uiMAX_REG_Y)
  \param uiNrBins Number of greybins for histogram ("dynamic range")
  \param fCliplimit Normalized cliplimit (higher values give more contrast)
  \return An error code, 0 meaning no error.
*/
template<typename T_Pixel>
inline int ClaheImpl<T_Pixel>::execute
(
  T_Pixel* pImage,
  unsigned int uiXRes,
  unsigned int uiYRes,
	T_Pixel Min,
  T_Pixel Max,
  unsigned int uiNrX,
  unsigned int uiNrY,
	unsigned int uiNrBins,
  float fCliplimit
)
{
  /* counters */
  unsigned int uiX, uiY;
  /* size of context. reg. and subimages */
  unsigned int uiXSize, uiYSize, uiSubX, uiSubY;
  /* auxiliary variables interpolation routine */
  unsigned int uiXL, uiXR, uiYU, uiYB;
  /* clip limit and region pixel count */
  unsigned long ulClipLimit, ulNrPixels;
  /* pointer to image */
  T_Pixel* pImPointer;
  
  /* lookup table used for scaling of input image */
  unsigned int uiNR_OF_GREY = _pow(2, (sizeof(T_Pixel) * 8));
  unsigned int* aLUT = new unsigned int[uiNR_OF_GREY];
  
  /* pointer to histogram and mappings*/
  unsigned long* pulHist, *pulMapArray;
  /* auxiliary pointers interpolation */
  unsigned long* pulLU, *pulLB, *pulRU, *pulRB;
  
  
  /* --------------------------------------------------------- */
  /* Just asserts parameter checking for speed in release mode */
  /* --------------------------------------------------------- */
  /* # of regions x-direction too large */
  assert(uiNrX <= uiMAX_REG_X);
  /* # of regions y-direction too large */
  assert(uiNrY <= uiMAX_REG_Y);
  /* x-resolution no multiple of uiNrX */
  assert(uiXRes % uiNrX == 0);
  /* y-resolution no multiple of uiNrY */
  assert(uiYRes % uiNrY == 0);
  /* maximum too large */
  assert(Max < uiNR_OF_GREY);
  /* minimum equal or larger than maximum */
  assert(Min < Max);
  /* at least 4 contextual regions required */
  assert(uiNrX >= 2 && uiNrY >= 2);
  /* is OK, immediately returns original image. */
  assert(fCliplimit != 1.0);
  /* default value when not specified */
  assert(uiNrBins > 0);
  /* Number of gray values must be a greater or equal to the number of bins. */
  assert((1 + (Max - Min)) >= uiNrBins);
  
  
#if 0
  /* # of regions x-direction too large */
  if (uiNrX > uiMAX_REG_X) return -1;
  /* # of regions y-direction too large */
  if (uiNrY > uiMAX_REG_Y) return -2;
  /* x-resolution no multiple of uiNrX */
  if (uiXRes % uiNrX) return -3;
  /* y-resolution no multiple of uiNrY */
  if (uiYRes & uiNrY) return -4;
  /* maximum too large */
  if (Max >= uiNR_OF_GREY) return -5;
  /* minimum equal or larger than maximum */
  if (Min >= Max) return -6;
  /* at least 4 contextual regions required */
  if (uiNrX < 2 || uiNrY < 2) return -7;
  /* is OK, immediately returns original image. */
  if (fCliplimit == 1.0) return 0;
  /* default value when not specified */
  if (uiNrBins == 0) uiNrBins = 128;
#endif
  
  
  pulMapArray = (unsigned long*)
                     malloc(sizeof(unsigned long) * uiNrX * uiNrY * uiNrBins);
  /* Not enough memory! (try reducing uiNrBins) */
  if (pulMapArray == 0) return -8;
  /* Initialize all the histograms at once */
  ::memset(pulMapArray, 0, sizeof(unsigned long) * uiNrX * uiNrY * uiNrBins);
  
  /* Actual size of contextual regions */
  uiXSize = uiXRes / uiNrX;
  uiYSize = uiYRes / uiNrY;
  ulNrPixels = (unsigned long)uiXSize * (unsigned long)uiYSize;
  
  /* Calculate actual cliplimit	 */
  if(fCliplimit > 0.0)
  {
    ulClipLimit = (unsigned long) (fCliplimit * (uiXSize * uiYSize) / uiNrBins);
    ulClipLimit = (ulClipLimit < 1UL) ? 1UL : ulClipLimit;
  }
  else
  {
    /* Large value, do not clip (AHE) */
    ulClipLimit = 1UL<<14;
  }
  
  
  /* Make lookup table for mapping of greyvalues */
  _makeLut(aLUT, Min, Max, uiNrBins);
  
  /* Calculate greylevel mappings for each contextual region */
  for (uiY = 0, pImPointer = pImage; uiY < uiNrY; uiY++)
  {
    for (uiX = 0; uiX < uiNrX; uiX++, pImPointer += uiXSize)
    {
      pulHist = &pulMapArray[uiNrBins * (uiY * uiNrX + uiX)];
      _makeHistogram(pImPointer, uiXRes, uiXSize, uiYSize,
                     pulHist, uiNrBins, aLUT);
      _clipHistogram(pulHist, uiNrBins, ulClipLimit);
      _mapHistogram(pulHist, Min, Max, uiNrBins, ulNrPixels);
    }
    /* skip lines, set pointer */
    pImPointer += (uiYSize - 1) * uiXRes;
  }
  
  
  /* Interpolate greylevel mappings to get CLAHE image */
  for (pImPointer = pImage, uiY = 0; uiY <= uiNrY; uiY++)
  {
    /* special case: top row */
    if (uiY == 0)
    {
      uiSubY = uiYSize >> 1;
      uiYU = 0;
      uiYB = 0;
    }
    else
    {
      /* special case: bottom row */
      if (uiY == uiNrY)
      {
        uiSubY = uiYSize >> 1;
        uiYU = uiNrY-1;
        uiYB = uiYU;
      }
      else
      {
        /* default values */
        uiSubY = uiYSize;
        uiYU = uiY - 1;
        uiYB = uiYU + 1;
      }
    }
    
    for (uiX = 0; uiX <= uiNrX; uiX++)
    {
      /* special case: left column */
      if (uiX == 0)
      {
        uiSubX = uiXSize >> 1;
        uiXL = 0;
        uiXR = 0;
      }
      else
      {
        /* special case: right column */
        if (uiX == uiNrX)
        {
          uiSubX = uiXSize >> 1;
          uiXL = uiNrX - 1;
          uiXR = uiXL;
        }
        else
        {
          /* default values */
          uiSubX = uiXSize; uiXL = uiX - 1; uiXR = uiXL + 1;
        }
      }
      
      pulLU = &pulMapArray[uiNrBins * (uiYU * uiNrX + uiXL)];
      pulRU = &pulMapArray[uiNrBins * (uiYU * uiNrX + uiXR)];
      pulLB = &pulMapArray[uiNrBins * (uiYB * uiNrX + uiXL)];
      pulRB = &pulMapArray[uiNrBins * (uiYB * uiNrX + uiXR)];
      _interpolate(pImPointer, uiXRes, pulLU, pulRU, pulLB, pulRB, uiSubX,
                  uiSubY, aLUT);
      
      /* set pointer on next matrix */
      pImPointer += uiSubX;
    }
    pImPointer += (uiSubY - 1) * uiXRes;
  }
  
  /* free space for histograms */
  free(pulMapArray);
  
  /* return status OK */
  return 0;
}


/**
  Compute the integral power of an integral base.
  FIXME: Use so math lib for this!
*/
template<typename T_Pixel>
inline unsigned int ClaheImpl<T_Pixel>::_pow
(
  unsigned int nBase,
  unsigned int nExp
)
{
  unsigned int nIdx = 1;
  unsigned int nRes = nBase;
  for ( ; nIdx < nExp; nIdx++)
  {
    nRes *= nBase;
  }
  return nRes;
}


/**
  FIXME: This will spin forever for very small values of the clipping limit,
         that is when there is a lot of pixels in excess, i.e. more to
         redistribute than available space ...
  
  Performs clipping of the histogram and redistribution of bins.
  The histogram is clipped and the number of excess pixels is counted.
  Afterwards the excess pixels are equally redistributed across the whole
  histogram (providing the bin count is smaller than the cliplimit).
  
  This method works by side-effect and will modify the contents of pulHistogram.
  
  \param pulHistogram The histogram. There will be \see{uiNrGreylevels} bins
    and so will be the useful length of this array. Every slot in the array
    represents the number of cases (pixels) for the bin.
  \param uiNrGreylevels The number of bins in the histogram.
  \param ulClipLimit The clip limit or, equivalently, the maximum slope of the
    cumulative histogram.
*/
template<typename T_Pixel>
inline void ClaheImpl<T_Pixel>::_clipHistogram
(
  unsigned long* pulHistogram,
  unsigned int uiNrGreylevels,
  unsigned long ulClipLimit
)
{
  unsigned long* pulBinPointer, *pulEndPointer, *pulHisto;
  unsigned long ulNrExcess, ulUpper, ulBinIncr, ulStepSize;
  long lBinExcess;
  
  /* calculate total number of excess pixels */
  ulNrExcess = 0; 
  pulBinPointer = pulHistogram;
  unsigned int nIdx = 0;
  for ( ; nIdx < uiNrGreylevels; nIdx++)
  {
    lBinExcess = (long) pulBinPointer[nIdx] - (long) ulClipLimit;
    if (lBinExcess > 0) ulNrExcess += lBinExcess;	  /* excess in current bin */
  }
  
  /* Second part: clip histogram and redistribute excess pixels in each bin */
  ulBinIncr = ulNrExcess / uiNrGreylevels;
      /* average binincrement */
      /* This will be our first approximation for redistribution */
  ulUpper =  ulClipLimit - ulBinIncr;
      /* Bins in between ulClipLimit and ulUpper will be set to cliplimit */
  
  nIdx = 0;
  for ( ; nIdx < uiNrGreylevels; nIdx++)
  {
    if (pulHistogram[nIdx] > ulClipLimit)
    {
      /* clip bin */
      pulHistogram[nIdx] = ulClipLimit;
    }
    else
    {
      if (pulHistogram[nIdx] > ulUpper)
      {
        /* high bin count */
        /* Everything in between clip limit and "upper" is given less than */
        /* than average so that we do not fill the bin to more than clip */
        /* limit which would not make sense. But we do give enough to reach */
        /* out the clip limit. */
	      ulNrExcess -= pulHistogram[nIdx] - ulUpper;
        pulHistogram[nIdx] = ulClipLimit;
      }
      else
      {
        /* low bin count */
        /* We give this guy the average */
	      ulNrExcess -= ulBinIncr;
        pulHistogram[nIdx] += ulBinIncr;
      }
    }
  }
  
  /* Finally, Redistribute remaining excess  */
  while (ulNrExcess)
  {
    pulEndPointer = &pulHistogram[uiNrGreylevels];
    pulHisto = pulHistogram;
    
    while (ulNrExcess && pulHisto < pulEndPointer)
    {
      /* If number in excess is less than the number of bins, we must skip */
      /* some bins. */
      ulStepSize = uiNrGreylevels / ulNrExcess;
	    if (ulStepSize < 1)
      {
        /* stepsize at least 1 */
        ulStepSize = 1;
      }
      
      /* Actual redistribution */
	    for ( pulBinPointer=pulHisto;
            pulBinPointer < pulEndPointer && ulNrExcess;
            pulBinPointer += ulStepSize)
      {
        if (*pulBinPointer < ulClipLimit)
        {
          (*pulBinPointer)++;	
          ulNrExcess--; /* reduce excess */
        }
      }
      
      /* Restart redistributing on other bin location. Yes but why ? */
	    pulHisto++;
    }
  }
  
}



/**
  Classifies the greylevels present in the array image into a greylevel histogram.
  
  The pLookupTable specifies the relationship between the greyvalue of the pixel
  (typically between 0 and 4095) and the corresponding bin in the histogram
  (usually containing only 128 bins).
  
  \param pImage The one channel 2D image as one flat array. The rows are layed
    out one after the other.
  \param uiXRes Complete image width.
  \param uiSizeX Image Region width
  \param uiSizeY Image Region height
  \param pulHistogram The histogram. There will be \see{uiNrGreylevels} bins
    and so will be the useful length of this array. Every slot in the array
    represents the number of cases (pixels) for the bin.
  \param uiNrGreylevels The number of bins in the histogram.
  \param pLookupTable The look up table. The index is a pixel input value and
           the value is the corresponding bin number in the histogram.
*/
template<typename T_Pixel>
inline void ClaheImpl<T_Pixel>::_makeHistogram
(
  T_Pixel* pImage,
  unsigned int uiXRes,
	unsigned int uiSizeX,
  unsigned int uiSizeY,
	unsigned long* pulHistogram,
	unsigned int uiNrGreylevels,
  unsigned int* pLookupTable
)
{
  T_Pixel* pImagePointer;
  unsigned int nIdx = 0;
  
  /* clear histogram */
  /* Not useful anymore. It is initialized from outside. */
  /* for (i = 0; i < uiNrGreylevels; i++) */
  /*  pulHistogram[i] = 0L; */
  
#if 0
  for ( ; nIdx < uiSizeY; nIdx++)
  {
    /* Points to the end of the current row for this region (one past the last one) */
    pImagePointer = &pImage[uiSizeX];
    
    /* Put all the pixels of this row for this region in the histogram */
    while (pImage < pImagePointer)
      pulHistogram[pLookupTable[*pImage++]]++;
    
    /* Skip to the beginning of the next line for this region  */
    /* No matter where the region is, we always have to go the complete image */
    /* width forward. Well not exactly since pImage is at the end of the */
    /* region where here. Therefore we need to go complete image width minus */
    /* the region size. */
    
    /* I think that Zuiderveld's implementation is not very clear. */
    /* That's why I substituted it for mine. */
    /*
    pImagePointer += uiXRes;
    pImage = &pImagePointer[-uiSizeX];
    */
    pImage += (uiXRes - uiSizeX);
  }
#endif
  
  /*
    Optimization :
    
    -1- The nIdx counter is useless. We just have to compute the last pixel
        to process once at the beginning and to compare with that. So we
        don't need to increment a counter each time we loop over the lines.
    -2- Compute uiXRes - uiSizeX once and for all outside the loop.
    -3- First loop comparison _must_ be useless.
    -4- Last loop image incrementation is useless. Better do loop condition
        before.
    
    So the loop becomes:
  */
  T_Pixel* pImageEnd = pImage + ((uiSizeY-1) * uiXRes) + uiSizeX;
  const unsigned int uiInc = uiXRes - uiSizeX;
  for ( ; ; )
  {
    pImagePointer = &pImage[uiSizeX];
    while (pImage < pImagePointer)
      pulHistogram[pLookupTable[*pImage++]]++;
    if (pImage == pImageEnd)
      break;
    pImage += uiInc;
  }
}



/**
  Calculates the equalized lookup table (mapping) by
  cumulating the input histogram. Note: lookup table is rescaled in
  range [Min..Max].
  
  This in effect computes the cumulative histogram so that the histogram
  itself will play the role of the mapping function afterwards.
  
  \param pulHistogram The histogram. There will be \see{uiNrGreylevels} bins
    and so will be the useful length of this array. Every slot in the array
    represents the number of cases (pixels) for the bin.
  \param Min The desired minimum gray value
  \param Max The desired maximum gray value
  \param uiNrGreylevels The number of bins in the histogram.
  \param ulNrOfPixels The complete number of pixels in the image
*/
template<typename T_Pixel>
inline void ClaheImpl<T_Pixel>::_mapHistogram
(
  unsigned long* pulHistogram,
  T_Pixel Min,
  T_Pixel Max,
  unsigned int uiNrGreylevels,
  unsigned long ulNrOfPixels
)
{
#if 0
  unsigned int i; 
  unsigned long ulSum = 0;
  const float fScale = ((float)(Max - Min)) / ulNrOfPixels;
  const unsigned long ulMin = (unsigned long) Min;
  
  for (i = 0; i < uiNrGreylevels; i++)
  {
    ulSum += pulHistogram[i];
    pulHistogram[i] = (unsigned long) (ulMin+ulSum*fScale);
    if (pulHistogram[i] > Max)
      pulHistogram[i] = Max;
  }
#endif

  /*
    Optimization:
    
    -1- Once the Max is reached, we don't need to compute anything anymore.
        Just set all the rest to Max
    
  */
  unsigned int i; 
  unsigned long ulSum = 0;
  const float fScale = ((float)(Max - Min)) / ulNrOfPixels;
  for (i = 0; i < uiNrGreylevels; i++)
  {
    ulSum += pulHistogram[i];
    pulHistogram[i] = (unsigned long) (Min+ulSum*fScale);
    
    if (pulHistogram[i] > Max)
    {
      /* Just set everything else to Max and break the loop */
      pulHistogram[i] = Max;
      i++;
      for ( ; i < uiNrGreylevels; i++)
      {
        pulHistogram[i] = Max;
      }
      break;
    }
  }
}



/**
  To speed up histogram clipping, the input image [Min,Max] is scaled down to
  [0,uiNrBins-1]. This function calculates the LUT (Look Up Table).
  
  \param pLUT The allocated LUT we need to fill up.
  \param nMin Desired minimum grey value.
  \param nMax Desired maximum grey value.
  \param uiNrBins Desired number of bins.
*/
template<typename T_Pixel>
inline void ClaheImpl<T_Pixel>::_makeLut
(
  unsigned int* pLUT,
  T_Pixel nMin,
  T_Pixel nMax,
  unsigned int uiNrBins
)
{
#if 1

  const T_Pixel BinSize = (T_Pixel) (1 + (nMax - nMin) / uiNrBins);
  int i;
  for (i = nMin; i <= nMax; i++)
    pLUT[i] = (i - nMin) / BinSize;

#else
  
  
  /*
    Optimization:
    
    -1- Division is just a performance pain. Actually we do _NOT_ have to
        compute anything. First bins numbers are 0, then 1, then 2, etc.
        We just need to know how many sequential input pixels will be assign
        to the same bin, which of course is just BinSize.
    
    -2- We do not need to check for the end of the pLut array since :
        
        uiNR_OF_GREY >= 1 + (Max - Min)                          i.e.
        uiNR_OF_GREY >= uiNrBins * (1 + (Max - Min)) / uiNrBins  i.e.
        uiNR_OF_GREY >= uiNrBins * BinSize                       
        
        That's good. But (1 + (Max - Min)) DOES NEED to be a multiple of
        uiNrBins otherwise some LUT slot won't be assigned any bin number.
        
        So, if we want to avoid this extra condition, we need to add to
        BinSize, in the case (1 + (Max - Min)) is not a multiple of uiNrBins,
        the necessary quantity during which the division in the old code will
        round up.
        
        That quantity of course is uiNrBins / 2. This is the
        limit starting from which we go to the next bin.
        
        Therefore this is a very good candidate to add to BinSize to get the
        actual real BinSize.
        
        But if we don't check for pLUT end, which we won't, that means that
        we may fill too much in the pLUT array (up to nNbGrayValues % uiNrBins)
        So we'll make an outer loop that will provide for the first uiNrBins - 1
        bins without checking for pLUT array end. Then for the last bin no,
        we'll check the pLut array end instead of filling for all the bin size.
    
    So the loop becomes:
  */
  const T_Pixel nNbGrayValues = (T_Pixel) (1 + (nMax - nMin));
  const T_Pixel BinSize = (T_Pixel) (nNbGrayValues / uiNrBins +
                                     uiNrBins / 2);
  unsigned int* pLUTEnd = pLUT + (nMax + 1);
  unsigned int nBinNo = 0;
  
  pLUT += nMin;
  for ( ; nBinNo < uiNrBins - 1; nBinNo++)
  {
    unsigned int* pLUTBinEnd = pLUT + BinSize;
    for ( ; pLUT < pLUTBinEnd; pLUT++)
      *pLUT = nBinNo;
  }
  
  /* Fill to the end for last bin no. */
  for ( ; pLUT < pLUTEnd; pLUT++)
    *pLUT = nBinNo;

#endif
}



/**
  Calculates the new greylevel assignments of pixels within a
  submatrix of the image with size uiXSize and uiYSize. This is done by a
  bilinear interpolation between four different mappings in order to eliminate
  boundary artifacts. It uses a division; since division is often an expensive
  operation, I added code to perform a logical shift instead when feasible.
  
  \param pImage pointer to input/output image; It's the beginning of the region
          to compute.
  \param uiXRes Complete width of the image
  \param pulMapLU Mapping function for the upper left region
  \param pulMapRU Mapping function for the upper right region
  \param pulMapLB Mapping function for the bottom left region
  \param pulMapRB Mapping function for the bottom right region
  \param uiXSize Width of the region
  \param uiYSize Height of the region
  \param pLUT lookup table containing mapping greyvalues to bins
*/
template<typename T_Pixel>
inline void ClaheImpl<T_Pixel>::_interpolate
(
  T_Pixel* pImage,
  unsigned int uiXRes,
  unsigned long* pulMapLU,
  unsigned long* pulMapRU,
  unsigned long* pulMapLB,
  unsigned long* pulMapRB,
  unsigned int uiXSize,
  unsigned int uiYSize,
  unsigned int* pLUT
)
{
  /* Pointer increment after processing row */
  const unsigned int uiIncr = uiXRes - uiXSize;
  T_Pixel GreyValue;
  
  /* Normalization factor */
  unsigned int uiNum = uiXSize * uiYSize;
  
  unsigned int uiXCoef = 0;
  unsigned int uiYCoef = 0;
  unsigned int uiXInvCoef = 0;
  unsigned int uiYInvCoef = 0;
  unsigned int uiShift = 0;
  
  /* If uiNum is not a power of two, use division */
  if (uiNum & (uiNum - 1))
  {
    uiYCoef = 0;
    uiYInvCoef = uiYSize;
    for
    (
      ;
      uiYCoef < uiYSize;
      uiYCoef++, uiYInvCoef--, pImage += uiIncr
    )
    {
      for
      (
        uiXCoef = 0, uiXInvCoef = uiXSize;
        uiXCoef < uiXSize;
        uiXCoef++, uiXInvCoef--
      )
      {
        /* get histogram bin value */
        GreyValue = pLUT[*pImage];
        
        *pImage++ = (T_Pixel)
        
        (
          
          (
             uiYInvCoef *
             (uiXInvCoef*pulMapLU[GreyValue] + uiXCoef * pulMapRU[GreyValue])
           +
             uiYCoef *
             (uiXInvCoef * pulMapLB[GreyValue] + uiXCoef * pulMapRB[GreyValue])
          )
          
          / uiNum
        );
      }
    }
  }
  else
  {
    /* avoid the division and use a right shift instead */
    while (uiNum >>= 1) uiShift++; /* Calculate 2log of uiNum */
    
    for
    (
      uiYCoef = 0, uiYInvCoef = uiYSize;
      uiYCoef < uiYSize;
	    uiYCoef++, uiYInvCoef--, pImage+=uiIncr
    )
    {
      for
      (
        uiXCoef = 0, uiXInvCoef = uiXSize;
        uiXCoef < uiXSize;
        uiXCoef++, uiXInvCoef--
      )
      {
        /* get histogram bin value */
        GreyValue = pLUT[*pImage];
        
        *pImage++ = (T_Pixel)
        
        (
          (
            uiYInvCoef *
            (
              uiXInvCoef * pulMapLU[GreyValue]
              +
              uiXCoef * pulMapRU[GreyValue]
            )
            
            +
            
            uiYCoef *
            (
              uiXInvCoef * pulMapLB[GreyValue]
              +
              uiXCoef * pulMapRB[GreyValue]
            )
          )
          >> uiShift
        );
        
      }
    }
  }
}


/**
  Default ctor.
*/
template<typename T_ItkImage>
inline ClaheITK<T_ItkImage>::ClaheITK()
  : _nMin(0),
    _nMax( (T_Pixel) -1),
    _uiNrX(8),
    _uiNrY(8),
    _uiNrBins(256),
    _fCliplimit(10.0),
    _bAutoAdaptToNbRegions(true),
    _pInput(NULL),
    _pFlatImage(NULL),
    _oClaheImpl(),
    _pImportFilter()
{
  _pImportFilter = T_ImportFilter::New();
}


/**
  Dtor.
*/
template<typename T_ItkImage>
inline ClaheITK<T_ItkImage>::~ClaheITK()
{
  free(_pFlatImage);
}


/**
  Set input image source.
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::SetInput(T_ItkImagePointer pInput)
{
  this->_pInput = pInput;
}


/**
  Get the output image.
*/
template<typename T_ItkImage>
inline typename T_ItkImage::Pointer ClaheITK<T_ItkImage>::GetOutput()
{
  assert(_pInput.GetPointer() != NULL);
  return _pImportFilter->GetOutput();
}


/**
  Update the image processing.
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::Update()
{
  assert(_pInput.GetPointer() != NULL);
  this->execute(_pInput);
  //_pImportFilter->Update();
}


/**
  Executes Clahe image processing.
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::execute(T_ItkImagePointer pImage)
{
  //
  // -1- Translate from ITK to flat array
  //
  const T_RegionType& roLargestRegion = pImage->GetLargestPossibleRegion();
  const T_ItkImage::SizeType& roLargestSize = roLargestRegion.GetSize();
  unsigned int uiXRes = (unsigned int) roLargestSize[0];
  unsigned int uiYRes = (unsigned int) roLargestSize[1];
  unsigned int nNbPixels = uiXRes * uiYRes;
  unsigned int nNbBytes = sizeof(T_Pixel) * nNbPixels;
  
  if (_pFlatImage != NULL)
  {
    free(_pFlatImage);
  }
  _pFlatImage = _fromItkToFlatArray(pImage, roLargestRegion, nNbBytes);
  
  //
  // -2- Apply CLAHE upon the flat array
  //
  T_ItkImagePointer pOutputImage = NULL;
  int nErrCode = 0;
  if (_bAutoAdaptToNbRegions)
  {
    nErrCode = _oClaheImpl.wrapedExecute
    (
      _pFlatImage, uiXRes, uiYRes, _nMin, _nMax, _uiNrX, _uiNrY,
      _uiNrBins, _fCliplimit
    );
  }
  else
  {
    nErrCode = _oClaheImpl.execute
    (
      _pFlatImage, uiXRes, uiYRes, _nMin, _nMax, _uiNrX, _uiNrY,
      _uiNrBins, _fCliplimit
    );
  }
  
  
  //
  // -3- Translate back from flat array to CLAHE
  //
  if (nErrCode == 0)
  {
    _fromFlatArrayToItk(_pFlatImage, roLargestSize, nNbPixels);
  }
}


/**
  Set the output minimum gray tone value.
  
  Output gray values shall be rescaled in nMin:nMax
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::setGrayLevelMin(T_Pixel nMin)
{
  _nMin = nMin;
}


/**
  Set the output maximum gray tone value.
  
  Output gray values shall be rescaled in nMin:nMax
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::setGrayLevelMax(T_Pixel nMax)
{
  _nMax = nMax;
}


/**
  Set the number of regions along which to divide the image along the X direction.
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::setNbRegionsX(unsigned int uiNrX)
{
  _uiNrX = uiNrX;
}


/**
  Set the number of regions along which to divide the image along the Y direction.
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::setNbRegionsY(unsigned int uiNrY)
{
  _uiNrY = uiNrY;
}


/**
  Set the number of bins in the histogram.
  
  Should of course be less than the number of pixels in a region.
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::setNbBins(unsigned int uiNrBins)
{
  _uiNrBins = uiNrBins;
}


/**
  Set the clip limit of the histogram or equivalently the slope of the
  cumulative histogram.
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::setCliplimit(float fCliplimit)
{
  _fCliplimit = fCliplimit;
}


/**
  Set whether we should adapt the image to the number of regions or not.
  
  If not, the width must be an integral multiple of the number of regions
  along the X axis and the height must be an integral multiple of the
  number of regions along the Y axis.
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::setAutoAdaptToNbRegions
(
  bool bAutoAdaptToNbRegions
)
{
  _bAutoAdaptToNbRegions = bAutoAdaptToNbRegions;
}


/**
  Get the output minimum gray tone value.
  
  Output gray values shall be rescaled in nMin:nMax
*/
template<typename T_ItkImage>
inline typename ClaheITK<T_ItkImage>::T_Pixel
ClaheITK<T_ItkImage>::getGrayLevelMin() const
{
  return _nMin;
}


/**
  Get the output maximum gray tone value.
  
  Output gray values shall be rescaled in nMin:nMax
*/
template<typename T_ItkImage>
inline typename ClaheITK<T_ItkImage>::T_Pixel
ClaheITK<T_ItkImage>::getGrayLevelMax() const
{
  return _nMax;
}


/**
  Get the number of regions along which to divide the image along the X direction.
*/
template<typename T_ItkImage>
inline unsigned int ClaheITK<T_ItkImage>::getNbRegionsX() const
{
  return _uiNrX;
}


/**
  Get the number of regions along which to divide the image along the Y direction.
*/
template<typename T_ItkImage>
inline unsigned int ClaheITK<T_ItkImage>::getNbRegionsY() const
{
  return _uiNrY;
}


/**
  Get the number of bins in the histogram.
  
  Should of course be less than the number of pixels in a region.
*/
template<typename T_ItkImage>
inline unsigned int ClaheITK<T_ItkImage>::getNbBins() const
{
  return _uiNrBins;
}


/**
  Get the clip limit of the histogram or equivalently the slope of the
  cumulative histogram.
*/
template<typename T_ItkImage>
inline float ClaheITK<T_ItkImage>::getCliplimit() const
{
  return _fCliplimit;
}


/**
  Get whether we should adapt the image to the number of regions or not.
  
  If not, the width must be an integral multiple of the number of regions
  along the X axis and the height must be an integral multiple of the
  number of regions along the Y axis.
*/
template<typename T_ItkImage>
inline bool ClaheITK<T_ItkImage>::getAutoAdaptToNbRegions() const
{
  return _bAutoAdaptToNbRegions;
}


/**
  Translate back from flat array to ITK.
  
  \param pFlatImage The gray tone image as a flat array.
  \param roImgSize The image size
  \param nNbPixels Complete number of pixels in the image.
  \return A smart pointer to the newly created 
*/
template<typename T_ItkImage>
inline void ClaheITK<T_ItkImage>::_fromFlatArrayToItk
(
  T_Pixel* pFlatImage,
  const T_SizeType& roImgSize,
  unsigned int nNbPixels
)
{
  T_ImportFilter::IndexType start;
  start.Fill(0);
  T_RegionType oRegionForImportFilter;
  oRegionForImportFilter.SetIndex(start);
  oRegionForImportFilter.SetSize(roImgSize);
  _pImportFilter->SetRegion(oRegionForImportFilter);
  
  double origin[2] = { 0.0, 0.0 };
  double spacing[2] = { 1.0, 1.0 };
  _pImportFilter->SetOrigin(origin);
  _pImportFilter->SetSpacing(spacing);
  const bool bImportImageFilterWillOwnTheBuffer = false; // Our shit to delete.
  _pImportFilter->SetImportPointer(pFlatImage, nNbPixels,
                                   bImportImageFilterWillOwnTheBuffer);
}


/**
  Translate ITK image to flat array.
  
  \return The flat C array constructed on the heap with malloc. It will be the
          responsibility of the caller to eventually deallocate this memory.
*/
template<typename T_ItkImage>
inline typename ClaheITK<T_ItkImage>::T_Pixel*
ClaheITK<T_ItkImage>::_fromItkToFlatArray
(
  T_ItkImagePointer pImage,
  const T_RegionType& roRegion,
  unsigned int nNbBytes
)
{
  typedef typename itk::ImageRegionConstIterator<T_ItkImage> T_ImgKIt;
  
  T_Pixel* pFlatImage = (T_Pixel*) malloc(nNbBytes);
  T_Pixel* itFlatArray = pFlatImage;
  
  T_ImgKIt itInput(pImage, roRegion);
  for (itInput.GoToBegin() ; !itInput.IsAtEnd(); ++itInput, itFlatArray++)
  {
    *itFlatArray = itInput.Get();
  }
  
  return pFlatImage;
}


// #define __clahe_h__
#endif

