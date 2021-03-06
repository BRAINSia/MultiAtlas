/*=========================================================================
  Copyright (c) IDEA LAB, UNC-Chapel Hill, 2013.

     MABMIS (Multi-Atlas-Based Multi-Image Segmentation)

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
// for math
#include <vcl_iostream.h>
#include <vnl/vnl_random.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matlab_print.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_svd_economy.h>
#include <string>
#include <vector>

// basic itk
#include "itkVector.h"

// registration
#include "itkImageRegistrationMethod.h"
#include "itkSymmetricForcesDemonsRegistrationFilter.h"

// interpolator
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

// reader / writer
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// filter
#include "itkResampleImageFilter.h"

#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkWarpVectorImageFilter.h"
#include "itkInverseDisplacementFieldImageFilter.h"

// for affine transformation
#include "itkTransform.h"
#include "itkAffineTransform.h"
#include "itkImageRegistrationMethod.h"

// for Diffeomorphic Demons
#include <itkCommand.h>
#include <itkDiffeomorphicDemonsRegistrationFilter.h>
#include <itkMultiResolutionPDEDeformableRegistration.h>

// including itksys::SystemTools::MakeDirectory(char*)
#include <itksys/SystemTools.hxx>
#include <metaCommand.h>

// To include all related header files
#include "IGR3D_MABMIS_TestingCLP.h"
#include "itkMABMISImageOperationFilter.h"
#include "itkMABMISDisplacementFieldFilter.h"
#include "itkMABMISSimulateData.h"
#include "itkMABMISImageRegistrationFilter.h"
#include "itkMABMISTreeOperation.h"
#include "itkMABMISBasicOperationFilter.h"

#include "itkMABMISAtlasXMLFile.h"
#include <algorithm>
static std::string ReplacePathSepForOS( const std::string & input )
{
  std::string output = input;
#ifdef _WIN32
  std::replace(output.begin(), output.end(), '/', FILESEP);
#else
  std::replace(output.begin(), output.end(), '\\', FILESEP);
#endif
  return output;
}

static std::string
GetExtension(const std::string & filename)
{
  std::string fileExt( itksys::SystemTools::GetFilenameLastExtension(filename) );
  //If the last extension is .gz, then need to pull off 2 extensions.
  //.gz is the only valid compression extension.
  if ( fileExt == std::string(".gz") )
    {
    fileExt = itksys::SystemTools::GetFilenameLastExtension( 
              itksys::SystemTools::GetFilenameWithoutLastExtension(filename) );
    fileExt += ".gz";
    }
  return ( fileExt );
}

static std::string
GetRootName(const std::string & filename)
{
  const std::string fileExt = GetExtension(filename);

  // Create a base filename
  // i.e Image.hdr --> Image
  if ( fileExt.length() > 0                    //Ensure that an extension was
                                               // found
       && filename.length() > fileExt.length() //Ensure that the filename does
                                               // not contain only the extension
       )
    {
    const std::string::size_type it = filename.find_last_of(fileExt);
    const std::string            baseName( filename, 0, it - ( fileExt.length() - 1 ) );
    return ( baseName );
    }
  //Default to return same as input when the extension is nothing (Analyze)
  return ( filename );
}

typedef double CoordinateRepType;

// basic data type
typedef unsigned char                                  CharPixelType;  // for image IO usage
typedef float                                          FloatPixelType; // for
typedef int                                            IntPixelType;
typedef short                                          ShortPixelType;
typedef float                                          InternalPixelType; // for internal processing usage
typedef itk::Vector<InternalPixelType, ImageDimension> VectorPixelType;

// basic image type
typedef itk::Image<CharPixelType, ImageDimension>     CharImageType;
typedef itk::Image<IntPixelType, ImageDimension>      IntImageType;
typedef itk::Image<ShortPixelType, ImageDimension>    ShortImageType;
typedef itk::Image<FloatPixelType, ImageDimension>    FloatImageType;
typedef itk::Image<InternalPixelType, ImageDimension> InternalImageType;
typedef itk::Image<VectorPixelType, ImageDimension>   DisplacementFieldType;

// basic iterator type
typedef itk::ImageRegionIterator<DisplacementFieldType> DisplacementFieldIteratorType;
typedef itk::ImageRegionIterator<InternalImageType>    InternalImageIteratorType;
typedef itk::ImageRegionIterator<CharImageType>        CharImageIteratorType;

// basic image reader/writer related type
typedef itk::ImageFileReader<CharImageType>     CharImageReaderType;
typedef itk::ImageFileReader<InternalImageType> InternalImageReaderType;
typedef itk::ImageFileWriter<InternalImageType> InternalImageWriterType;

typedef itk::WarpImageFilter<InternalImageType, InternalImageType, DisplacementFieldType> InternalWarpFilterType;
typedef itk::ImageFileWriter<CharImageType>                                              CharImageWriterType;
typedef itk::ImageFileWriter<IntImageType>                                               IntImageWriterType;
typedef itk::ImageFileWriter<FloatImageType>                                             FloatImageWriterType;
typedef itk::ImageFileWriter<ShortImageType>                                             ShortImageWriterType;

typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFieldReaderType;
typedef itk::ImageFileWriter<DisplacementFieldType> DisplacementFieldWriterType;

//////////////////////////////////////////////////////////////////////////////
// image filter type
typedef itk::ResampleImageFilter<InternalImageType, InternalImageType>          ResampleFilterType;
typedef itk::HistogramMatchingImageFilter<InternalImageType, InternalImageType> InternalHistMatchFilterType;

////////////////////////////////////////////////////////////////////////////
// operation on deformation fields
typedef itk::WarpVectorImageFilter<DisplacementFieldType, DisplacementFieldType,
                                   DisplacementFieldType>                                      WarpVectorFilterType;
typedef itk::InverseDisplacementFieldImageFilter<DisplacementFieldType,
                                                DisplacementFieldType>
  InverseDisplacementFieldImageFilterType;
typedef itk::AddImageFilter<DisplacementFieldType, DisplacementFieldType, DisplacementFieldType> AddImageFilterType;

// global bool variables to adjust the  procedure

bool isEvaluate = false; // if false, we do not know the ground-truth of labels
bool isDebug = false;     // false;//true; // if true, print out more information

int localPatchSize = 1; // (2r+1)*(2r+1)*(2r+1) is the volume of local patch

// demons registration parameters
// int iterInResolutions[4][3]={{5,3,2},{10,5,5},{15,10,5},{20,15,10}};
// int itereach = 2; //
// int itereach0 = 0;int itereach1 = 1;int itereach2 = 2;int itereach3 = 3;
// double sigmaDef = 1.5;
// double sigmaDef10 = 1.0;double sigmaDef15 = 1.5;double sigmaDef20 = 2.0;
// double sigmaDef25 = 2.5;double sigmaDef30 = 3.0;double sigmaDef35 = 3.5;
bool doHistMatch = true;
bool doHistMatchTrue = true; bool doHistMatchFalse = false;

typedef itk::Statistics::MABMISSimulateData<InternalImageType, InternalImageType> DataSimulatorType;
DataSimulatorType::Pointer datasimulator;

typedef itk::Statistics::MABMISImageOperationFilter<CharImageType, CharImageType> ImageOperationFilterType;
ImageOperationFilterType::Pointer imgoperator;
typedef itk::Statistics::MABMISDisplacementFieldFilter<InternalImageType,
                                                      InternalImageType> DisplacementFieldOperationFilterType;
DisplacementFieldOperationFilterType::Pointer dfoperator;
typedef itk::Statistics::MABMISImageRegistrationFilter<CharImageType, CharImageType> ImageRegistrationFilterType;
ImageRegistrationFilterType::Pointer regoperator;
typedef itk::Statistics::MABMISTreeOperation<InternalImageType, InternalImageType> TreeOperationType;
TreeOperationType::Pointer treeoperator;

typedef itk::Statistics::MABMISBasicOperationFilter<CharImageType, CharImageType> BasicOperationFilterType;
BasicOperationFilterType::Pointer basicoperator;

//

typedef itk::Vector<ShortPixelType, ImageDimension>      ShortVectorPixelType;
typedef itk::Image<ShortVectorPixelType, ImageDimension> ShortDisplacementFieldType;
typedef itk::ImageFileWriter<ShortDisplacementFieldType>  ShortDisplacementFieldWriterType;

void ReadImgInfo(std::vector<std::string> imgfilenames);

void HistogramMatching(InternalImageType::Pointer inputImage, InternalImageType::Pointer referenceImage,
                       InternalImageType::Pointer & outputImage);

void TreeBasedRegistrationFastOniTree(vnl_vector<int> itree,          // the incremental tree
                                      int root,                       // the tree root
                                      int itree_size,                 // the tree size
                                      itk::MABMISAtlas* atlasTree,    // the atlas tree
                                      itk::MABMISImageData*imageData, // the images to be segmented
                                      std::vector<int> iterations,    // registration parameters
                                      double sigma);

void LabelFusion(std::string curSampleImgName, std::string outSampleSegName, std::vector<std::string>  allWarpedAtlasImgNames,
                 std::vector<std::string> allDisplacementFieldNames, std::vector<std::string> allAtlasSegNames, int numOfAtlases);

void GaussianWeightedLabelFusion(InternalImageType::Pointer curSampleImgPtr, InternalImageType::Pointer outSampleSegPtr,
                                 InternalImageType::Pointer* warpedImgPtrs, InternalImageType::Pointer* warpedSegPtrs,
                                 int numOfAtlases);

void RegistrationOntoTreeRoot(vnl_vector<int> itree,          // the incremental tree
                              int root,                       // the tree root
                              int itree_size,                 // the tree size
                              itk::MABMISAtlas* atlasTree,    // the atlas tree
                              itk::MABMISImageData*imageData, // the images to be segmented
                              std::vector<int> iterations,    // registration parameters
                              double sigma);

// void PairwiseRegistrationOnTreeViaRoot(int root, int filenumber, int atlas_size, std::vector<std::string> sub_ids);
void PairwiseRegistrationOnTreeViaRoot(int root, itk::MABMISImageData* imageData, itk::MABMISAtlas* atlasTree,
                                       std::vector<int> iterations, double sigma);

// void MultiAtlasBasedSegmentation(int filenumber, int atlas_size, std::vector<std::string> sub_ids);
int MultiAtlasBasedSegmentation(int root, itk::MABMISImageData* imageData, itk::MABMISAtlas* atlasTree);

int Testing(itk::MABMISImageData* imageData, itk::MABMISAtlas* atlasTree,
            std::vector<int> iterations, double sigma)
{
  // Validate the tree size is correct
  if( atlasTree->m_TreeSize != atlasTree->m_Tree.size() ||
      atlasTree->m_TreeSize != atlasTree->m_NumberAllAtlases )
    {
    std::cerr << "ERROR: The atlas tree size is incorrect!!" << std::endl;
    return -1;
    }
  /////////////////////////////////////
  // the atlas tree
  int tree_size = atlasTree->m_NumberAllAtlases;

  // incremental tree: the atlas tree PLUS images to be segmented
  int             itree_size = atlasTree->m_NumberAllAtlases + imageData->m_NumberImageData;
  vnl_vector<int> itree(itree_size);
  for( int i = 0; i < itree_size; ++i )
    {
    itree.put(i, 0);
    }
  // copy current tree into an incremental tree
  for( int i = 0; i < tree_size; ++i )
    {
    itree.put(i, atlasTree->m_Tree[i]);                                   // combinative tree: all atlases including
                                                                          // simulated data
    }
  treeoperator->FindRoot(itree, atlasTree->m_TreeSize, atlasTree->m_TreeRoot);
  treeoperator->GetTreeHeight(itree, atlasTree->m_TreeSize, atlasTree->m_TreeHeight);
  int root = atlasTree->m_TreeRoot;  // root node
  // Validate the tree is correct
  for( int n = 0; n < atlasTree->m_TreeSize; ++n )
    {
    bool valid = true;
    if( atlasTree->m_Tree[n] >= atlasTree->m_TreeSize )
      {
      valid = false;
      }
    if( atlasTree->m_Tree[n] < 0 && n != atlasTree->m_TreeRoot )
      {
      valid = false;
      }
    if( !valid )
      {
      std::cerr << "ERROR: The atlas tree size is incorrect!!" << std::endl;
      return -1;
      }
    }

  if( isDebug )
    {
    treeoperator->SaveTreeWithInfo(itree, itree_size, imageData->m_DataDirectory + "itree.txt");
    }

  // build the incremental tree
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Build incremental tree ... " << std::endl;

  int            totalNumAtlases = atlasTree->m_NumberAllAtlases;
  int            totalNumFiles = totalNumAtlases + imageData->m_NumberImageData;
  std::vector<std::string> allfilenames(totalNumFiles);
  for( int i = 0; i < totalNumFiles; ++i )
    {
		
    if( i < totalNumAtlases )
      {
      allfilenames[i] = ReplacePathSepForOS(atlasTree->m_AtlasDirectory + atlasTree->m_AtlasFilenames[i]);
      }
    else
      {
      allfilenames[i] = ReplacePathSepForOS(imageData->m_DataDirectory + imageData->m_ImageFileNames[i - totalNumAtlases]);
      }
    }

  std::cout << "1. Calculating pairwise distance cross training and test images ..." << std::endl;
  vnl_matrix<double> cross_dist(imageData->m_NumberImageData, totalNumFiles);
  for( int i = 0; i < imageData->m_NumberImageData; ++i )
    {
    for( int j = 0; j < totalNumFiles; ++j )
      {
      cross_dist.put(i, j, imgoperator->calculateDistanceMSD(allfilenames[totalNumAtlases + i], allfilenames[j]) );
      }
    }

  if( isDebug )
    {
    basicoperator->SaveMatrix2File(cross_dist, imageData->m_NumberImageData, totalNumFiles,
                                   imageData->m_DataDirectory + "sample2atlas_dist.txt");
    }

  // indices to for atlas files and image files to be segmented.
  int* atlas_cur = new int[totalNumFiles];
  for( int i = 0; i < totalNumFiles; ++i )
    {
    atlas_cur[i] = 0;
    }
  for( int i = 0; i < totalNumAtlases; ++i )
    {
    atlas_cur[i] = i;
    }
  int* images_cur = new int[imageData->m_NumberImageData];
  for( int i = 0; i < imageData->m_NumberImageData; ++i )
    {
    images_cur[i] = i + totalNumAtlases;
    }

  // build the incremental tree
  treeoperator->SetAllAtlasNumber(totalNumAtlases);
  std::cout << "2. Build the tree... " << std::endl;
  itree = treeoperator->BuildIncrementalTree(imageData->m_NumberImageData, // #images to be segmented
                                             totalNumAtlases,              // #atlases, including original ones and
                                                                           // simulated ones
                                             totalNumFiles,                // total #image files, =
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           //
                                                                           // imageData->m_NumberImageData+totalNumAtlases
                                             images_cur,                   // file indices for images to be segmented
                                             atlas_cur,                    // file indices for atlas images
                                             cross_dist,                   // pairwise distance between images
                                             itree);                       // the output incremental tree

  // delete
  delete[] atlas_cur;
  delete[] images_cur;

  // check if new tree is actually a tree
  const bool iTreeValid = treeoperator->ValidateTree(itree, itree_size);
  if( !iTreeValid )
    {
    std::cout << "The new built tree is NOT valid!" << std::endl;

    return 105; // EXIT_FAILURE;
    }
  // save new tree
  if( isDebug )
    {
    treeoperator->SaveTreeWithInfo(itree, itree_size, imageData->m_DataDirectory + "itree.txt");
    // std::cout << "Tree is built up!" << std::endl;
    }
  std::cout << "Done! " << std::endl;

  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Register all images to the root... " << std::endl;
  // register all images to the root, if not done yet.
  RegistrationOntoTreeRoot(itree,      // the incremental tree
                           root,       // the tree root
                           itree_size, // the tree size
                           atlasTree,
                           imageData,
                           iterations, sigma // registration parameters
                           );
  std::cout << "Done. " << std::endl;

  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Pairwise registration via the tree root... " << std::endl;
  PairwiseRegistrationOnTreeViaRoot(root, imageData, atlasTree, iterations, sigma);
  std::cout << "Done. " << std::endl;

  ////////////////////////////////
  // do all multi-atlas based segmentation with weighted label fusion
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "Segment images using label fusion... " << std::endl;
  int numIter = MultiAtlasBasedSegmentation(root, imageData, atlasTree);

  std::cout << "Done. " << std::endl;

  // do all multi-atlas based segmentation with weighted label fusion
  /////////////////////////////////
  //
  // basicoperator->RemoveFile("*_*_cbq_000*");
  // basicoperator->RemoveFile("*_*_deform_000*");
  //
  /////////////////////////////////

  // remove intermediate results
  // deformedImageFileName and deformation fields
  if( !isDebug )
    {
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "Removing intermediate registration files..." << std::endl;
    for( int n = 0; n < imageData->m_NumberImageData; ++n )
      {
      char n_str[10], m_str[10];
      sprintf(n_str, "%03d", n);
      for( int m = 0; m < atlasTree->m_NumberAllAtlases - atlasTree->m_NumberSimulatedAtlases; ++m )
        {
        sprintf(m_str, "%03d", m);

        std::string imHdr = std::string(m_str) + "_to_" + "testImage" + n_str + "_cbq_000.nii.gz";
        std::string dfMha = std::string(m_str) + "_to_" + "testImage" + n_str + "_deform_000.nii.gz";
        imHdr = imageData->m_DataDirectory + imHdr;
        dfMha = imageData->m_DataDirectory + dfMha;
        basicoperator->RemoveFile(imHdr.c_str() );
        basicoperator->RemoveFile(dfMha.c_str() );

        if( m == root )
          {
          imHdr = std::string("testImage") + n_str + "_to_" + m_str + "_cbq_reg.nii.gz";
          dfMha = std::string("testImage") + n_str + "_to_" + m_str + "_deform_000.nii.gz";
          imHdr = imageData->m_DataDirectory + imHdr;
          dfMha = imageData->m_DataDirectory + dfMha;
          basicoperator->RemoveFile(imHdr.c_str() );
          basicoperator->RemoveFile(dfMha.c_str() );
          }
        }
      for( int m = 0; m < imageData->m_NumberImageData; ++m )
        {
        if( m == n )
          {
          continue;
          }
        sprintf(m_str, "%03d", m);
        std::string imHdr = std::string("testImage") + n_str + "_to_testImage" + m_str + "_cbq_000.nii.gz";
        std::string dfMha = std::string("testImage") + n_str + "_to_testImage" + m_str + "_deform_000.nii.gz";
        imHdr = imageData->m_DataDirectory + imHdr;
        dfMha = imageData->m_DataDirectory + dfMha;
        basicoperator->RemoveFile(imHdr.c_str() );
        basicoperator->RemoveFile(dfMha.c_str() );
        }
      }
    }
  // remove segmentation results in iterations.
  for( int n = 0; n < imageData->m_NumberImageData; ++n )
    {
    const std::string imageFileName = imageData->m_DataDirectory + imageData->m_ImageFileNames[n];
    const std::string baseFileName = GetRootName(imageFileName);
    char        i_str[10];
    bool        copied = false;
    for( int i = numIter - 1; i >= 0; i-- )
      {
      sprintf(i_str, "%03d", i);
      const std::string segHdr = baseFileName + "_seg_" + i_str + ".nii.gz";
      if( itksys::SystemTools::FileExists(segHdr.c_str(), true) )
        {
        if( !copied )
          {
          std::string segHdr_save = baseFileName + "_seg.nii.gz";

          const size_t dir_sep = segHdr_save.find_last_of(FILESEP);
          if( dir_sep != std::string::npos )
            {
            segHdr_save = imageData->m_OutputDirectory + segHdr_save.substr(dir_sep + 1, std::string::npos);
            }
          else
            {
            segHdr_save = imageData->m_OutputDirectory + segHdr_save;
            }

          itksys::SystemTools::CopyFileAlways(segHdr.c_str(), segHdr_save.c_str() );
          copied = true;
          }
        basicoperator->RemoveFile(segHdr.c_str() );
        }
      }
    }

  std::cout << "Done. " << std::endl;
  return EXIT_SUCCESS;
}

template <class T>
int DoIt( itk::MABMISImageData* imageData, itk::MABMISAtlas* atlasTree,
          std::vector<int> iterations, double sigma)
{
  return Testing(imageData, atlasTree, iterations, sigma);
}

int main( int argc, char *argv[] )
{
  PARSE_ARGS;

  // step 1: read the test image list
  itk::MABMISImageDataXMLFileReader::Pointer imageListXMLReader = itk::MABMISImageDataXMLFileReader::New();
  ImageListXML = ReplacePathSepForOS(ImageListXML); 
  AtlaseTreeXML = ReplacePathSepForOS(AtlaseTreeXML); 
  OutputFolder = ReplacePathSepForOS(OutputFolder); 

  imageListXMLReader->SetFilename(ImageListXML);

  try
    {
    imageListXMLReader->GenerateOutputInformation();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading file" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  itk::MABMISImageData * inputImageData = imageListXMLReader->GetOutputObject();

  // if the data path is empty, use the path of the xml file instead
  if( inputImageData->m_DataDirectory.size() <= 1 )
    {
    const size_t dir_sep = ImageListXML.find_last_of(FILESEP);
    if( dir_sep != std::string::npos )
      {
      inputImageData->m_DataDirectory = ImageListXML.substr(0, dir_sep);
      }
    else
      {
      inputImageData->m_DataDirectory.resize(0);
      }
    }
  if( !inputImageData->m_DataDirectory.empty() )
    {
    if( !(inputImageData->m_DataDirectory[inputImageData->m_DataDirectory.size() - 1] == FILESEP) )
      {
      inputImageData->m_DataDirectory = inputImageData->m_DataDirectory + FILESEP;
      }
    }

  // output directory
  if( OutputFolder.empty() )
    {
    inputImageData->m_OutputDirectory = inputImageData->m_DataDirectory;
    }
  else
    {
    if( !(OutputFolder[OutputFolder.size() - 1] == FILESEP) )
      {
      OutputFolder = OutputFolder + FILESEP;
      }
    inputImageData->m_OutputDirectory = OutputFolder;
    }

  // load the tree-structured atlases that is generated in the training step
  itk::MABMISAtlasXMLFileReader::Pointer treeAtlasXMLReader = itk::MABMISAtlasXMLFileReader::New();
  treeAtlasXMLReader->SetFilename(AtlaseTreeXML);
  try
    {
    treeAtlasXMLReader->GenerateOutputInformation();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading file" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  itk::MABMISAtlas * atlasTree = treeAtlasXMLReader->GetOutputObject();

  // set the atlas path as the same path as the xml file.
  const size_t dir_sep = AtlaseTreeXML.find_last_of(FILESEP);
  if( dir_sep != std::string::npos )
    {
    atlasTree->m_AtlasDirectory = ReplacePathSepForOS(AtlaseTreeXML.substr(0, dir_sep + 1) + atlasTree->m_AtlasDirectory);
    }

  if( atlasTree->m_AtlasDirectory.size() == 0 )
    {
    atlasTree->m_AtlasDirectory = ReplacePathSepForOS(".");
    }

  if( !(atlasTree->m_AtlasDirectory[atlasTree->m_AtlasDirectory.size() - 1] == FILESEP) )
    {
    atlasTree->m_AtlasDirectory = ReplacePathSepForOS(atlasTree->m_AtlasDirectory + FILESEP );
    }

  // Will look into getting rid of these global variables later ---Xiaofeng
  datasimulator = DataSimulatorType::New();
  imgoperator = ImageOperationFilterType::New();
  dfoperator = DisplacementFieldOperationFilterType::New();
  regoperator = ImageRegistrationFilterType::New();
  treeoperator = TreeOperationType::New();
  basicoperator = BasicOperationFilterType::New();

  int retVal = DoIt<unsigned short>( inputImageData, atlasTree, iterations, SmoothingKernelSize);

  delete inputImageData;
  delete atlasTree;
  return retVal;
}

void LabelFusion(std::string curSampleImgName, std::string outSampleSegName, std::vector<std::string> allWarpedAtlasImgNames,
                 std::vector<std::string> allDisplacementFieldNames, std::vector<std::string> allAtlasSegNames, int numOfAtlases)
{
  // create pointers for each images including: cur_image, our_seg,
  // warpedImg(after histogram matching), warpedSeg (after warping)

  // cur sample image pointer
  InternalImageType::Pointer curSampleImgPtr = 0;

  imgoperator->ReadImage(curSampleImgName, curSampleImgPtr);

  // output sample label image pointer
  InternalImageType::Pointer outSampleSegPtr = 0;
  imgoperator->ReadImage(curSampleImgName, outSampleSegPtr);

  // atlas related images pointers
  InternalImageType::Pointer*  warpedImgPtrs = new InternalImageType::Pointer[numOfAtlases];
  InternalImageType::Pointer*  warpedSegPtrs = new InternalImageType::Pointer[numOfAtlases];
  InternalImageType::IndexType index;
  index[0] = 19; index[1] = 19; index[2] = 8;
  for( int i = 0; i < numOfAtlases; ++i )
    {
    // histogram match each warped atlas to the current sample
    warpedImgPtrs[i] = 0;
    InternalImageType::Pointer curWarpedAtlasPtr = 0;
    imgoperator->ReadImage(allWarpedAtlasImgNames[i], curWarpedAtlasPtr);
    HistogramMatching(curWarpedAtlasPtr, curSampleImgPtr, warpedImgPtrs[i]);

    // load each deformation field
    DisplacementFieldType::Pointer deformFieldPtr = 0;
    dfoperator->ReadDisplacementField(allDisplacementFieldNames[i], deformFieldPtr);

    // warp each atlas label image to the current sample space
    InternalImageType::Pointer curAtlasSegPtr = 0;
    imgoperator->ReadImage(allAtlasSegNames[i], curAtlasSegPtr);
    dfoperator->ApplyDisplacementField(curAtlasSegPtr, deformFieldPtr, warpedSegPtrs[i], false);
    }

  // weighted label fusion
  GaussianWeightedLabelFusion(curSampleImgPtr, outSampleSegPtr, warpedImgPtrs, warpedSegPtrs, numOfAtlases);

  // output segmentation image
  imgoperator->WriteImage(outSampleSegName, outSampleSegPtr);

  delete[] warpedImgPtrs;
  delete[] warpedSegPtrs;

  return;
}

void GaussianWeightedLabelFusion(InternalImageType::Pointer curSampleImgPtr, InternalImageType::Pointer outSampleSegPtr,
                                 InternalImageType::Pointer* warpedImgPtrs, InternalImageType::Pointer* warpedSegPtrs,
                                 int numOfAtlases)
{
  // introduce the iterators of each image
  // InternalImageIteratorType sampleImgIt(curSampleImgPtr,
  // curSampleImgPtr->GetLargestPossibleRegion());//GetRequestedRegion());
  InternalImageIteratorType sampleSegIt(outSampleSegPtr, outSampleSegPtr->GetLargestPossibleRegion() ); //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        //
                                                                                                        // GetRequestedRegion());

  InternalImageIteratorType* warpedSegIts = new InternalImageIteratorType[numOfAtlases];

  typedef itk::ConstNeighborhoodIterator<InternalImageType> InternalImageNeighborhoodIteratorType;
  InternalImageNeighborhoodIteratorType* warpedImgNeighborhoodIts =
    new InternalImageNeighborhoodIteratorType[numOfAtlases];
  InternalImageType::SizeType radius;
  radius[0] = localPatchSize; radius[1] = localPatchSize; radius[2] = localPatchSize;

  InternalImageNeighborhoodIteratorType sampleImgNeighborhoodIt(radius, curSampleImgPtr,
                                                                curSampleImgPtr->GetRequestedRegion() );

  InternalImageNeighborhoodIteratorType sampleImgNeighborIt(radius, curSampleImgPtr,
                                                            curSampleImgPtr->GetLargestPossibleRegion() );
  for( int i = 0; i < numOfAtlases; ++i )
    {
    // InternalImageIteratorType it1( warpedImgPtrs[i],
    // warpedImgPtrs[i]->GetLargestPossibleRegion());//GetRequestedRegion() );
    // warpedImgIts[i] = it1;

    InternalImageIteratorType it2( warpedSegPtrs[i], warpedSegPtrs[i]->GetLargestPossibleRegion() ); //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     // GetRequestedRegion()
                                                                                                     // );
    warpedSegIts[i] = it2;

    InternalImageNeighborhoodIteratorType neighborIt(radius, warpedImgPtrs[i],
                                                     warpedImgPtrs[i]->GetLargestPossibleRegion() );
    warpedImgNeighborhoodIts[i] = neighborIt;
    }

  int maxLabel = 0; // the maximum labels
  // find out the maximum label from the label images
  for( int  i = 0; i < numOfAtlases; ++i )
    {
    warpedSegIts[i].GoToBegin();
    }
  while( !warpedSegIts[0].IsAtEnd() )
    {
    for( int  i = 0; i < numOfAtlases; ++i )
      {
      if( warpedSegIts[i].Get() > maxLabel )
        {
        maxLabel = warpedSegIts[i].Get();
        }
      ++warpedSegIts[i];
      }
    }


  int     temp1, temp2;
  double  kernel_sigma = 1; // std of Gaussian approximation, set to the median of all distances
  int*    index = new int[numOfAtlases];
  double* sort1 = new double[numOfAtlases];
  int*    label_index = new int[maxLabel + 1];

  int*    label_pool = new int[numOfAtlases];
  double* mse = new double[numOfAtlases];
  double* weight_sum =  new double[maxLabel + 1];
  // for each voxel do label fusion,(vx, vy, vz) is the current location
  double kernel_sigma_square = kernel_sigma * kernel_sigma;
  for( int i = 0; i < numOfAtlases; ++i )
    {
    warpedImgNeighborhoodIts[i].GoToBegin();
    warpedSegIts[i].GoToBegin();
    }
  for( sampleImgNeighborIt.GoToBegin(), sampleSegIt.GoToBegin();
       !sampleImgNeighborIt.IsAtEnd();
       ++sampleImgNeighborIt, ++sampleSegIt )
    {
    double msec = 0.0;
    for( int i = 0; i < numOfAtlases; ++i )
      {
      label_pool[i] = warpedSegIts[i].Get();
      for( int n = 0; n < sampleImgNeighborIt.Size(); ++n )
        {
        temp1 = sampleImgNeighborIt.GetPixel(n);
        temp2 = warpedImgNeighborhoodIts[i].GetPixel(n);
        // add together differences, squared sum
        msec += (temp1 - temp2) * (temp1 - temp2);
        }
      mse[i] = msec;
      }
    // sort the mean square difference
    for( int i = 0; i < numOfAtlases; ++i )
      {
      index[i] = i;
      sort1[i] = mse[i];
      }
    basicoperator->bubbleSort(sort1, index, numOfAtlases);
    kernel_sigma_square = sort1[(int)( (numOfAtlases - 1) / 2)] + 0.0001;
    // weight each label
    for( int i = 0; i < maxLabel + 1; ++i )
      {
      weight_sum[i] = 0.0;
      }
    for( int i = 0; i < numOfAtlases; ++i )
      {
      // add-up all weights
      weight_sum[label_pool[i]] += exp(0.0 - (mse[i]) / (2.0 * kernel_sigma_square) );
      }
    // weighted label fusion
    for( int i = 0; i < maxLabel + 1; ++i )
      {
      label_index[i] = i;
      }

    // determine the label with the maximum weight
    // basicoperator->bubbleSort(weight_sum, label_index, maxLabel+1);
    // int label_final;
    // label_final = label_index[maxLabel];
    int    label_final = 0;
    double weight_final = 0.0;
    for( int i = 0; i < maxLabel + 1; ++i )
      {
      if( weight_sum[i] > weight_final )
        {
        label_final = i;
        weight_final = weight_sum[i];
        }
      }

    // assign label to the segmentations
    sampleSegIt.Set(label_final);
    for( int i = 0; i < numOfAtlases; ++i )
      {
      ++warpedImgNeighborhoodIts[i];
      ++warpedSegIts[i];
      }
    }
  delete[] label_pool;
  delete[] mse;
  delete[] weight_sum;
  //
  // delete[] warpedImgIts;
  delete[] warpedSegIts;
  delete[] warpedImgNeighborhoodIts;
  delete[] sort1;
  delete[] index;
  delete[] label_index;

  return;
}

void TreeBasedRegistrationFastOniTree(vnl_vector<int> itree,          // the incremental tree
                                      int root,                       // the tree root
                                      int itree_size,                 // the tree size
                                      itk::MABMISAtlas* atlasTree,    // the atlas tree
                                      itk::MABMISImageData*imageData, // the images to be segmented
                                      std::vector<int> iterations,
                                      double sigma
                                      )
// (vnl_vector<int> itree_t, int itree_size_t, int atlas_size_t, int simulate_size_t, int sample_size_t, std::vector<std::string>
// originalIntensityImageFileNames, std::vector<std::string> deformedImgFileNames, std::vector<std::string> deformationFieldFileNames)
{
  // int root = -1; //
  // treeoperator->FindRoot(itree_t, itree_size_t, root);

  // calculate the depth of each node
  int* node_depth = new int[itree_size];

  for( int i = 0; i < itree_size; ++i )
    {
    int curnode = i;
    node_depth[i] = 0;
    bool isRoot = false;

    while( !isRoot )
      {
      if( itree[curnode] >= 0 )
        {
        node_depth[i]++;
        curnode = itree[curnode];
        }
      else
        {
        isRoot = true;
        }
      }
    }

  // sort all nodes by its depth
  int* index = new int[itree_size];
    {
    //
    double* arr = new double[itree_size];
    for( int i = 0; i < itree_size; ++i )
      {
      arr[i] = node_depth[i];
      index[i] = i;
      }
    basicoperator->bubbleSort(arr, index, itree_size); // index[0] = root
    delete[] arr;
    }

  int atlas_image_size = atlasTree->m_NumberAllAtlases - atlasTree->m_NumberSimulatedAtlases;
  int atlas_total_size = atlasTree->m_NumberAllAtlases;

  // start to register each image to the root node step by step
  for( int ii = 1; ii < itree_size; ++ii ) // starting from 1, since index[0] showing the root
    {
    // if (isDebug)
    //	std::cout << ii << ", ";
    int curnode = index[ii];
    int parentnode = itree[curnode];

    std::cout << "Between nodes: [" << curnode << ", " << parentnode << "]" << std::endl;

    if( (curnode >= atlas_image_size) && (curnode < atlas_total_size) )
      {
      // this is a simulated image, then pass
      continue;
      }

    // first check whether registrations from atlas images to root exist. If exist, no registration is needed
    bool isRegistered = false;
    char i_str[10], root_str[10], p_str[10];
    sprintf(i_str, "%03d", curnode);
    sprintf(root_str, "%03d", root);
    sprintf(p_str, "%03d", parentnode);

    std::string dfFileName, df_ImageFileName;
    if( curnode == root )
      {
      isRegistered = true;
      }
    else if( curnode < atlas_image_size )
      {
      dfFileName = std::string(i_str) + "_to_" + root_str + "_deform_000.nii.gz";
      df_ImageFileName = std::string(i_str) + "_to_" + root_str + "_cbq_reg.nii.gz";
      dfFileName = atlasTree->m_AtlasDirectory + dfFileName;
      df_ImageFileName = atlasTree->m_AtlasDirectory + df_ImageFileName;

      if( (itksys::SystemTools::FileExists(dfFileName.c_str(), true) &&
           itksys::SystemTools::FileExists(df_ImageFileName.c_str(), true) )
          )
        {
        isRegistered = true;
        }
      }
    else if( curnode >= atlas_total_size )
      {
      int  n = curnode - atlas_total_size;
      char n_str[10]; sprintf(n_str, "%03d", n);
      dfFileName = std::string("testImage") + n_str + "_to_" + root_str + "_deform_000.nii.gz";
      df_ImageFileName = std::string("testImage") + n_str + "_to_" + root_str + "_cbq_reg.nii.gz";

      dfFileName = imageData->m_DataDirectory + dfFileName;
      df_ImageFileName = imageData->m_DataDirectory + df_ImageFileName;

      isRegistered = false;
      }
    // parent deformation fields
    std::string df_ParentFileName;
    if( parentnode < atlas_image_size )
      {
      df_ParentFileName = std::string(p_str) + "_to_" + root_str + "_deform_000.nii.gz";
      df_ParentFileName = atlasTree->m_AtlasDirectory + df_ParentFileName;
      }
    else if( parentnode >= atlas_total_size )
      {
      int  n = parentnode - atlas_total_size;
      char n_str[10]; sprintf(n_str, "%03d", n);
      df_ParentFileName = std::string("testImage") + n_str + "_to_" + root_str + "_deform_000.nii.gz";
      df_ParentFileName = imageData->m_DataDirectory + df_ParentFileName;
      }

    // do registration
    if( isRegistered )
      {
      continue; // goto next one
      }

    std::string rootImageFile = ReplacePathSepForOS(atlasTree->m_AtlasDirectory + atlasTree->m_AtlasFilenames[root]);
    std::string testImageFile;
    if( curnode < atlas_image_size )
      {
      testImageFile = ReplacePathSepForOS(atlasTree->m_AtlasDirectory + atlasTree->m_AtlasFilenames[curnode]);    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    // atlas
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    //
                                                                                                                    // image
      }
    else
      {
      testImageFile = imageData->m_DataDirectory + imageData->m_ImageFileNames[curnode - atlas_total_size];   // test
                                                                                                              // image
      }
    if( parentnode == root ) // direct registration without initial deformation
      {
      regoperator->DiffeoDemonsRegistrationWithParameters(
        rootImageFile,
        testImageFile,
        df_ImageFileName,
        dfFileName,
        sigma, doHistMatch, iterations);
      }
    else // registration with initial deformation
      {
      if( (parentnode < atlas_image_size) || (parentnode >= atlas_total_size) )
        {
        // parent is a real image
        regoperator->DiffeoDemonsRegistrationWithInitialWithParameters(
          rootImageFile,
          testImageFile,
          df_ParentFileName,
          df_ImageFileName,
          dfFileName,
          sigma, doHistMatch, iterations);
        }
      else
        {
        // parent is a simulated image

        std::string index_string;  basicoperator->myitoa( parentnode - atlas_image_size, index_string, 3 );
        std::string sim_df_FileName = "simulated_deform_" + index_string + ".nii.gz";
        sim_df_FileName = atlasTree->m_AtlasDirectory + sim_df_FileName;

        regoperator->DiffeoDemonsRegistrationWithInitialWithParameters(
          rootImageFile,
          testImageFile,
          sim_df_FileName,
          df_ImageFileName,
          dfFileName,
          sigma, doHistMatch, iterations);
        }
      }
    }
  std::cout << "Done." << std::endl;

  delete[] node_depth;
  delete[] index;
}

// other filter
// histogram matching
void HistogramMatching(InternalImageType::Pointer inputImage,
                       InternalImageType::Pointer referenceImage,
                       InternalImageType::Pointer & outputImage)
{
  InternalHistMatchFilterType::Pointer matcher = InternalHistMatchFilterType::New();

  matcher->SetInput( inputImage );
  matcher->SetReferenceImage( referenceImage );
  matcher->SetNumberOfHistogramLevels( 1024 ); // 1024
  matcher->SetNumberOfMatchPoints( 7 );
  matcher->ThresholdAtMeanIntensityOn();
  try
    {
    matcher->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << err << std::endl;
    return;
    }
  outputImage = matcher->GetOutput();
  return;
}

// void RegistrationOntoTreeRoot(vnl_vector<int> itree, int root, int filenumber, int itree_size, int atlas_size, int
// simulate_size, int sample_size, std::vector<std::string> sub_ids, std::vector<std::string>  segfilenames)
void RegistrationOntoTreeRoot(vnl_vector<int> itree,          // the incremental tree
                              int root,                       // the tree root
                              int itree_size,                 // the tree size
                              itk::MABMISAtlas* atlasTree,    // the atlas tree
                              itk::MABMISImageData*imageData, // the images to be segmented
                              std::vector<int> iterations,    // registration parameters
                              double sigma
                              )
{
  // do registration on the tree
  std::cout << "Start registration ... " << std::endl;

  // to get the numbers original defined in the data;
  const int atlas_image_size = atlasTree->m_NumberAllAtlases - atlasTree->m_NumberSimulatedAtlases;
  const int test_image_size = imageData->m_NumberImageData;

// do tree-based registration
  TreeBasedRegistrationFastOniTree(itree,            // the incremental tree
                                   root,             // the tree root
                                   itree_size,       // the tree size
                                   atlasTree,        // the atlas tree
                                   imageData,        // the images to be segmented
                                   iterations, sigma // registration parameters
                                   );

  // end of do registration to the root

  //
  // basicoperator->RemoveFile("simulated_*");

  std::cout << "Reverse deformation field ... " << std::endl;;
  // reverse each deformation field (from the root)
  std::string rootImageTag;
  char        root_str[10], i_str[10];
  sprintf(root_str, "%03d", root);
  rootImageTag = root_str;
  std::string fixedImageTag;

  for( int i = 0; i < atlas_image_size + test_image_size; i++ )
    {
    //if( isDebug )
      {
      std::cout << i << ", ";
      }
    if( i == root )
      {
      continue;
      }

    // prepare file names

    // std::string originalImgImageFileName = sub_ids[root] + "_cbq_000.nii.gz";
    // std::string originalSegImageFileName = sub_ids[root] + "_seg_000.nii.gz";
    std::string rootImageFileName = ReplacePathSepForOS(
        atlasTree->m_AtlasDirectory + atlasTree->m_AtlasFilenames[root] );
    std::string rootSegmentFileName = ReplacePathSepForOS(
        atlasTree->m_AtlasDirectory + atlasTree->m_AtlasSegmentationFilenames[root] );

    std::string fixedImageFileName;
    std::string deformedImageFileName;
    std::string deformedImageFileNameImg;
    std::string deformedSegmentFileName;

    std::string deformationFileName;
    std::string invDeformationFileName;

    if( i < atlas_image_size )
      {
      fixedImageFileName = ReplacePathSepForOS(atlasTree->m_AtlasDirectory + atlasTree->m_AtlasFilenames[i]);
      sprintf(i_str, "%03d", i);
      fixedImageTag = i_str;

      deformedImageFileName = atlasTree->m_AtlasDirectory + rootImageTag + "_to_" + fixedImageTag + +"_cbq_000.nii.gz";
      deformedImageFileNameImg = atlasTree->m_AtlasDirectory + rootImageTag + "_to_" + fixedImageTag + +"_cbq_000.nii.gz";

      deformedSegmentFileName = atlasTree->m_AtlasDirectory + rootImageTag + "_to_" + fixedImageTag + +"_seg_000.nii.gz";

      deformationFileName = atlasTree->m_AtlasDirectory + fixedImageTag + "_to_" + rootImageTag + "_deform_000.nii.gz";
      invDeformationFileName = atlasTree->m_AtlasDirectory + rootImageTag + "_to_" + fixedImageTag + "_deform_000.nii.gz";
      }
    else
      {
      fixedImageFileName = imageData->m_DataDirectory + imageData->m_ImageFileNames[i - atlas_image_size];
      sprintf(i_str, "%03d", i - atlas_image_size);
      fixedImageTag = std::string("testImage") + i_str;

      deformedImageFileName = imageData->m_DataDirectory + rootImageTag + "_to_" + fixedImageTag + +"_cbq_000.nii.gz";
      deformedImageFileNameImg = imageData->m_DataDirectory + rootImageTag + "_to_" + fixedImageTag + +"_cbq_000.nii.gz";

      deformedSegmentFileName = imageData->m_DataDirectory + rootImageTag + "_to_" + fixedImageTag + +"_seg_000.nii.gz";

      deformationFileName = imageData->m_DataDirectory + fixedImageTag + "_to_" + rootImageTag + "_deform_000.nii.gz";
      invDeformationFileName = imageData->m_DataDirectory + rootImageTag + "_to_" + fixedImageTag + "_deform_000.nii.gz";
      }

    // if exist??
    bool isReversedExist = true;
    if( !itksys::SystemTools::FileExists(deformedImageFileName.c_str(), true) )
      {
      isReversedExist = false;
      }
    if( !itksys::SystemTools::FileExists(deformedImageFileNameImg.c_str(), true) )
      {
      isReversedExist = false;
      }
    if( !itksys::SystemTools::FileExists(invDeformationFileName.c_str(), true) )
      {
      isReversedExist = false;
      }
    if( isReversedExist )
      {
      continue;
      }

    // reverse deformation field
    DisplacementFieldType::Pointer deformationField = 0;
    dfoperator->ReadDisplacementField(deformationFileName, deformationField);

    DisplacementFieldType::Pointer inversedDisplacementField = DisplacementFieldType::New();
    inversedDisplacementField->SetRegions(deformationField->GetLargestPossibleRegion() );
    inversedDisplacementField->SetSpacing(deformationField->GetSpacing() );
    inversedDisplacementField->SetDirection(deformationField->GetDirection() );
    inversedDisplacementField->SetOrigin(deformationField->GetOrigin() );
    inversedDisplacementField->Allocate();

    dfoperator->InverseDisplacementField3D(deformationField, inversedDisplacementField);
    // output reversed deformation field

    dfoperator->WriteDisplacementField(invDeformationFileName, inversedDisplacementField);
    // //update
    regoperator->DiffeoDemonsRegistrationWithInitialWithParameters(
      fixedImageFileName, rootImageFileName,
      invDeformationFileName,
      deformedImageFileName, invDeformationFileName,
      sigma, doHistMatch, iterations);
        
    // CompressDisplacementField2Short(inversedDisplacementFieldFileName);
    // apply deformation field on seg file
    if( isEvaluate )
      {
      dfoperator->ApplyDisplacementFieldAndWriteWithFileNames(
        rootSegmentFileName,
        invDeformationFileName,
        deformedSegmentFileName, false);
      }

    std::cout << "Done!" << std::endl;
    } // end of reverse each deformation field (from the root)
      // std::cout << "done." << std::endl;
}

void PairwiseRegistrationOnTreeViaRoot(int root,
                                       itk::MABMISImageData* imageData,
                                       itk::MABMISAtlas* atlasTree,
                                       std::vector<int> iterations,
                                       double sigma)
{
  /////////////////////////////////
  // with the help of root node, obtain the pairwise registration among all images
  std::cout << "Generate all pairwise registration ..." << std::endl;;

  int atlas_image_size = atlasTree->m_NumberAllAtlases - atlasTree->m_NumberSimulatedAtlases;
  int test_image_size = imageData->m_NumberImageData;

  // do for all real images, including both atlases and test images
  for( int all_index = 0; all_index < atlas_image_size + test_image_size; all_index++ )
    {
    if( isDebug )
      {
      std::cout << all_index << ": ";
      }
    if( all_index == root )
      {
      continue;
      }

    char a_str[10], s_str[10], root_str[10];
    sprintf(root_str, "%03d", root);
    sprintf(a_str, "%03d", all_index);

    std::string movingImageFileName;
    std::string movingSegmentFileName;
    // get the moving filenames
    std::string movingImageTag;
    if( all_index < atlas_image_size )
      {
      movingImageFileName = ReplacePathSepForOS(atlasTree->m_AtlasDirectory + atlasTree->m_AtlasFilenames[all_index]);
      movingSegmentFileName = atlasTree->m_AtlasDirectory + atlasTree->m_AtlasSegmentationFilenames[all_index];
      movingImageTag = std::string(a_str);
      }
    else
      {
      movingImageFileName = imageData->m_DataDirectory + imageData->m_ImageFileNames[all_index - atlas_image_size];
      movingSegmentFileName = movingImageFileName;
      const std::string movingSegmentBaseFileName = GetRootName(movingSegmentFileName);

      movingSegmentFileName = movingSegmentBaseFileName + "_seg_000.nii.gz";
      sprintf(a_str, "%03d", all_index - atlas_image_size);
      movingImageTag = std::string("testImage") + a_str;
      }
    // for (int sample_index = atlas_size; sample_index < filenumber; sample_index++)
    for( int sample_index = 0; sample_index < imageData->m_NumberImageData; sample_index++ )
      {
      if( isDebug )
        {
        std::cout << "[" << all_index << "," << sample_index << "], ";
        }

      if( sample_index == all_index - atlas_image_size )
        {
        continue;
        }

      std::string fixedImageFileName = imageData->m_DataDirectory + imageData->m_ImageFileNames[sample_index];

      sprintf(s_str, "%03d", sample_index);
      std::string fixedImageTag = std::string("testImage") + s_str;

      std::string deformedImageFileName = movingImageTag + "_to_" + fixedImageTag + "_cbq_000.nii.gz";
      std::string deformedImageFileNameImg = movingImageTag + "_to_" + fixedImageTag + "_cbq_000.nii.gz";
      std::string deformedSegmentFileName = movingImageTag + "_to_" + "testImage" + s_str + "_seg_000.nii.gz";
      std::string deformedSegmentFileNameImg = movingImageTag + "_to_" + "testImage" + s_str + "_seg_000.nii.gz";
      std::string dfFileName = movingImageTag + "_to_" + fixedImageTag + "_deform_000.nii.gz";

      deformedImageFileName = imageData->m_DataDirectory + deformedImageFileName;
      deformedImageFileNameImg = imageData->m_DataDirectory + deformedImageFileNameImg;
      dfFileName = imageData->m_DataDirectory + dfFileName;

      deformedSegmentFileName = imageData->m_DataDirectory + deformedSegmentFileName;
      deformedSegmentFileNameImg = imageData->m_DataDirectory + deformedSegmentFileNameImg;

      // std::string fixedImgImageFileName = sub_ids[sample_index] + "_cbq_000.nii.gz";
      // std::string movingImgImageFileName = sub_ids[all_index] + "_cbq_000.nii.gz";
      // std::string movingSegImageFileName = sub_ids[all_index] + "_seg_000.nii.gz";

      // std::string deformedImgImageFileName = sub_ids[sample_index] + "_" + sub_ids[all_index] + "_cbq_000.nii.gz";
      // std::string deformedSegImageFileName = sub_ids[sample_index] + "_" + sub_ids[all_index] + "_seg_000.nii.gz";
      // std::string deformedImgImageFileNameImg = sub_ids[sample_index] + "_" + sub_ids[all_index] + "_cbq_000.nii.gz";
      // std::string deformedSegImageFileNameImg = sub_ids[sample_index] + "_" + sub_ids[all_index] + "_seg_000.nii.gz";
      // std::string deformationFieldFileName = sub_ids[sample_index] + "_" + sub_ids[all_index] + "_deform_000.nii.gz";

      // if exist
      bool isComposeExist = true;
      if( !itksys::SystemTools::FileExists(deformedImageFileName.c_str(), true) )
        {
        isComposeExist = false;
        }
      if( !itksys::SystemTools::FileExists(deformedImageFileNameImg.c_str(), true) )
        {
        isComposeExist = false;
        }

      if( !itksys::SystemTools::FileExists(dfFileName.c_str(), true) )
        {
        isComposeExist = false;
        }
      if( isComposeExist )
        {
        continue;
        }

      // compose
      // std::string inputDisplacementFieldFileName = sub_ids[root] + "_" + sub_ids[all_index] + "_deform_000.nii.gz";
      // std::string secondDisplacementFieldFileName = sub_ids[sample_index] + "_" + sub_ids[root] + "_deform_000.nii.gz";
      std::string inputDisplacementFieldFileName;
      if( all_index < atlas_image_size )
        {
        inputDisplacementFieldFileName = atlasTree->m_AtlasDirectory + movingImageTag + "_to_" + root_str
          + "_deform_000.nii.gz";
        }
      else
        {
        inputDisplacementFieldFileName  = imageData->m_DataDirectory + movingImageTag + "_to_" + root_str
          + "_deform_000.nii.gz";
        }
      std::string secondDisplacementFieldFileName = imageData->m_DataDirectory + root_str + "_to_" + fixedImageTag
        + "_deform_000.nii.gz";
      // output composed deformation field
      dfoperator->ComposeDisplacementFieldsAndSave(inputDisplacementFieldFileName,
                                                  secondDisplacementFieldFileName,
                                                  dfFileName);

      // update: refine the deformation
      regoperator->DiffeoDemonsRegistrationWithInitialWithParameters(
        fixedImageFileName,
        movingImageFileName,
        dfFileName,
        deformedImageFileName,
        dfFileName,
        sigma, doHistMatch, iterations);

      if( isEvaluate )
        {
        dfoperator->ApplyDisplacementFieldAndWriteWithFileNames(
          movingSegmentFileName,
          dfFileName,
          deformedSegmentFileName, false);
        }
      } // end of sample_index
    }   // end of all_index

  std::cout << "Done!" << std::endl;
}

// void MultiAtlasBasedSegmentation(int filenumber, int atlas_size, std::vector<std::string> sub_ids)
int MultiAtlasBasedSegmentation(int,
                                itk::MABMISImageData* imageData,
                                itk::MABMISAtlas* atlasTree)
{
  std::cout << "Start to process segmentation..." << std::endl;
  int iter = 1;

  bool isConverged = false;  // the iterative update can stop

  int atlas_image_size = atlasTree->m_NumberAllAtlases - atlasTree->m_NumberSimulatedAtlases;

  int test_image_size = imageData->m_NumberImageData;

  while( !isConverged )
    {
    // if (isDebug)
    std::cout << "iteration " << iter << ": ";
    // iter == 1, only using atlases; otherwise using all images
    // number of other image to be used
    int numOtherAtlases = 0;
    if( iter == 1 )
      {
      numOtherAtlases = atlas_image_size;
      }
    else
      {
      numOtherAtlases = atlas_image_size + test_image_size - 1;
      }
    // do label fusion for each new sample
    for( int sample_index = 0; sample_index < test_image_size; sample_index++ )
      {
      // if (isDebug)
      //	std::cout << sample_index << ":" << sub_ids[sample_index] << ", ";
      std::cout << sample_index << ", ";

      char s_str[10];
      sprintf(s_str, "%03d", sample_index);
      std::string testImageTag = std::string("testImage") + s_str;

      // the sample to be processed

      std::string testImageFileName = imageData->m_DataDirectory + imageData->m_ImageFileNames[sample_index];

      // std::string curSampleImgName = sub_ids[sample_index] +"_cbq_000.nii.gz";
      // the output label image of current sample
      std::string index_string;
      basicoperator->myitoa( iter, index_string, 3 );
      std::string testSementFileName = testImageFileName;
      const std::string testSegmentBaseFileName = GetRootName(testSementFileName);
      testSementFileName = testSegmentBaseFileName + "_seg_" + index_string + ".nii.gz";

      // std::string   outSampleSegName = sub_ids[sample_index] + "_seg_" + index_string + ".nii.gz";

      // prepare all current atlases names
      std::vector<std::string> allWarpedAtlasImgNames;
      std::vector<std::string> allDisplacementFieldNames;
      std::vector<std::string> allAtlasSegNames;

      int atlas_index = 0;
      for( int i = 0; i < atlas_image_size + test_image_size; ++i )
        {
        if( (iter == 1) && (i >= atlas_image_size) )
          {
          continue;
          }
        if( i - atlas_image_size == sample_index )
          {
          continue;
          }

        std::string movingImageTag;
        char        i_str[10];

        if( i < atlas_image_size )
          {
          sprintf(i_str, "%03d", i);
          movingImageTag = std::string(i_str);
          }
        else
          {
          sprintf(i_str, "%03d", i - atlas_image_size);
          movingImageTag = std::string("testImage") + i_str;
          }

        // assign names
        // std::string allWarpedAtlasImgName = sub_ids[sample_index] + "_" + sub_ids[i] + "_cbq_000.nii.gz";
        // std::string allDisplacementFieldName = sub_ids[sample_index] + "_" + sub_ids[i] + "_deform_000.nii.gz";
        std::string warpedImageFileName = movingImageTag + "_to_" + testImageTag + "_cbq_000.nii.gz";
        std::string dfFileName = movingImageTag + "_to_" + testImageTag + "_deform_000.nii.gz";

        warpedImageFileName = imageData->m_DataDirectory + warpedImageFileName;
        dfFileName = imageData->m_DataDirectory + dfFileName;

        std::string atlasSegName;
        if( i < atlas_image_size )
          {
          atlasSegName = atlasTree->m_AtlasDirectory + atlasTree->m_AtlasSegmentationFilenames[i];
          }
        else
          {
          atlasSegName = imageData->m_DataDirectory + imageData->m_ImageFileNames[i - atlas_image_size];
          const std::string atlasSegBaseName = GetRootName( atlasSegName );

          basicoperator->myitoa( iter - 1, index_string, 3 );
          atlasSegName = atlasSegBaseName + "_seg_" + index_string + ".nii.gz";
          }

        allWarpedAtlasImgNames.push_back(warpedImageFileName);
        allDisplacementFieldNames.push_back(dfFileName);
        allAtlasSegNames.push_back(atlasSegName);

        // go to next atlas
        atlas_index++;
        }
      // do weighted label fusion

      LabelFusion(testImageFileName,  // the input image to be segmented
                  testSementFileName, // the output segmented image
                  allWarpedAtlasImgNames,
                  allDisplacementFieldNames,
                  allAtlasSegNames,
                  numOtherAtlases);
      }
    std::cout << std::endl;

    // to decide if converged
    if( iter >= 3 )
      {
      isConverged = true; //
      }
    // calculate overlap rate
    if( isEvaluate )
      {
      // pairwise segmentation accuracy

      // pairwise segmentation consistency
      }
    iter++;
    }

  return iter;
}

