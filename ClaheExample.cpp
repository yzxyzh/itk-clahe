/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: BilateralImageFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2009-03-16 21:52:47 $
  Version:   $Revision: 1.24 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4996 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif


#include "clahe.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc != 9 ) 
  { 
    std::cerr << "Usage: " << std::endl
              << argv[0] << " inputImageFile outputImageFile GrayMin GrayMax "
              << std::endl
              << "NbRegionsX NbRegionsY NbBins ClipLimit"
              << std::endl;
    
    return EXIT_FAILURE;
  }
  
  typedef unsigned char InputPixelType;
  typedef unsigned char OutputPixelType;
  
  typedef itk::Image<InputPixelType, 2> InputImageType;
  typedef itk::Image<OutputPixelType, 2> OutputImageType;
  typedef itk::ImageFileReader<InputImageType>  ReaderType;
  typedef ClaheImageFilter<InputImageType> FilterType;
  
  FilterType* pFilter = new FilterType();
  
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  
  pFilter->setGrayLevelMin( static_cast<InputPixelType>(atoi(argv[3])) );
  pFilter->setGrayLevelMax( static_cast<InputPixelType>(atoi(argv[4])) );
  pFilter->setNbRegionsX( static_cast<unsigned int>(atoi(argv[5])) );
  pFilter->setNbRegionsY( static_cast<unsigned int>(atoi(argv[6])) );
  pFilter->setNbBins( static_cast<unsigned int>(atoi(argv[7])) );
  pFilter->setCliplimit( atof(argv[8]) );
  
  //
  // This filter beaks the streaming pipeline. Therefore we must "update"
  // the reader (i.e. execute it completely) before giving it to the reader.
  //
  reader->Update();
  pFilter->SetInput(reader->GetOutput());
  
  //
  // ... and we must "update" the filter (i.e. execute it completely) before
  // giving it to someone else, here the writer.
  //
  pFilter->Update();
  
  //
  // Write the output file.
  //
  typedef unsigned char WritePixelType;
  typedef itk::Image<WritePixelType, 2>        WriteImageType;
  typedef itk::ImageFileWriter< WriteImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( pFilter->GetOutput() );
  writer->Update();
  
  return EXIT_SUCCESS;
}
