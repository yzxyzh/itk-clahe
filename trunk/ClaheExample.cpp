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
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile outputImageFile GrayMin, GrayMax, NbRegionsX, NbRegionsY, NbBins, ClipLimit" << std::endl;
    return EXIT_FAILURE;
  }
  
  
  //  Software Guide : BeginLatex
  //
  //  The image types are instantiated using pixel type and dimension.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef    unsigned char    InputPixelType;
  typedef    unsigned char    OutputPixelType;

  typedef itk::Image< InputPixelType,  2 >   InputImageType;
  typedef itk::Image< OutputPixelType, 2 >   OutputImageType;
  // Software Guide : EndCodeSnippet


  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  
  
  typedef ClaheITK<InputImageType> FilterType;
  FilterType* filter = new FilterType();
  
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  
  filter->setGrayLevelMin( static_cast<InputPixelType>(atoi(argv[3])) );
  filter->setGrayLevelMax( static_cast<InputPixelType>(atoi(argv[4])) );
  filter->setNbRegionsX( static_cast<unsigned int>(atoi(argv[5])) );
  filter->setNbRegionsY( static_cast<unsigned int>(atoi(argv[6])) );
  filter->setNbBins( static_cast<unsigned int>(atoi(argv[7])) );
  filter->setCliplimit( atof(argv[8]) );
  
  //
  // This filter beaks the streaming pipeline. Therefore we must "update"
  // the reader before giving it to 
  reader->Update();
  filter->SetInput(reader->GetOutput());
  
  typedef unsigned char                          WritePixelType;
  typedef itk::Image< WritePixelType, 2 >        WriteImageType;
  typedef itk::ImageFileWriter< WriteImageType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  std::cerr << "================ argv[2] = " << argv[2] << std::endl;
  
  writer->SetFileName( argv[2] );
  
  filter->Update();
  writer->SetInput( filter->GetOutput() );
  writer->Update();
  
  return EXIT_SUCCESS;
}
