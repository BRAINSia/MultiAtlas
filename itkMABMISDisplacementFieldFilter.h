#ifndef __itkMABMISDisplacementFieldFilter_h
#define __itkMABMISDisplacementFieldFilter_h

#include <itkImage.h>
#include <itkImageToImageFilter.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>

#include "itkImage.h"
#include "itkWarpImageFilter.h"
#include "itkAddImageFilter.h"

// interpolator
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkMABMISImageOperationFilter.h"
#include "itkWarpVectorImageFilter.h"

#define ImageDimension 3

namespace itk
{
namespace Statistics
{
template <class TInputImage, class TOutputImage>
class MABMISDisplacementFieldFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef MABMISDisplacementFieldFilter                  Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  typedef TInputImage                 ImageType;
  typedef typename ImageType::Pointer ImagePointerType;

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

  typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFieldReaderType;
  typedef itk::ImageFileWriter<DisplacementFieldType> DisplacementFieldWriterType;

  typedef itk::WarpVectorImageFilter<DisplacementFieldType, DisplacementFieldType,
                                     DisplacementFieldType>         WarpVectorFilterType;
  typedef itk::AddImageFilter<DisplacementFieldType, DisplacementFieldType,
                              DisplacementFieldType>                AddImageFilterType;

// basic iterator type
  typedef itk::ImageRegionIterator<DisplacementFieldType> DisplacementFieldIteratorType;
  typedef itk::ImageRegionIterator<InternalImageType>    InternalImageIteratorType;
  typedef itk::ImageRegionIterator<CharImageType>        CharImageIteratorType;

// interpolator type
  typedef itk::LinearInterpolateImageFunction<InternalImageType, double>          InternalLinearInterpolatorType;
  typedef itk::NearestNeighborInterpolateImageFunction<InternalImageType, double> InternalNNInterpolatorType;

  typedef itk::WarpImageFilter<InternalImageType, InternalImageType, DisplacementFieldType> InternalWarpFilterType;

  typedef itk::Statistics::MABMISImageOperationFilter<ImageType, ImageType> ImageOperationType;
  typename ImageOperationType::Pointer imgoperator;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MABMISDisplacementFieldFilter, ImageToImageFilter);

  DisplacementFieldType::SpacingType   df_spacing;
  DisplacementFieldType::DirectionType df_direction;
  DisplacementFieldType::PointType     df_origin;

  int ReadDisplacementField(std::string filename, DisplacementFieldType::Pointer & deformationfield);

  void WriteDisplacementField(std::string  filename, DisplacementFieldType::Pointer deformationfield);

  void ComposeDisplacementFieldsAndSave(std::string inputDisplacementFieldFileName, std::string deformationFieldFileName,
                                       std::string composedDisplacementFieldFileName);

  void ComposeDisplacementFields(DisplacementFieldType::Pointer input, DisplacementFieldType::Pointer deformationField,
                                DisplacementFieldType::Pointer & composedDisplacementField);

  void InverseDisplacementField3D(DisplacementFieldType::Pointer deformationField,
                                 DisplacementFieldType::Pointer & deformationFieldInverse);

  void ApplyDisplacementField(InternalImageType::Pointer movingImage, DisplacementFieldType::Pointer deformationField,
                             InternalImageType::Pointer & deformedImage, bool isLinearInterpolator);

  void ApplyDisplacementFieldAndWriteWithFileNames(std::string movingImageName, std::string deformationFieldFileName,
                                                  std::string deformedImageName, bool isLinearInterpolator);

  void ApplyDisplacementFieldAndWriteWithTypeWithFileNames(std::string  movingImageFileName,
                                                          std::string deformationFieldFileName,
                                                          std::string deformedImageFileName, bool isLinear);

  void DownResampleDisplacementField(std::string deformationFieldFileName, std::string resampledDisplacementFieldFileName,
                                    int sampleRate);

  void UpResampleDisplacementField(std::string deformationFieldFileName, std::string  resampledDisplacementFieldFileName,
                                  int sampleRate);

  itkSetMacro(Imx, int);
  itkSetMacro(Imy, int);
  itkSetMacro(Imz, int);
private:
  MABMISDisplacementFieldFilter(const Self &); // purposely not implemented
  void operator=(const Self &);               // purposely not implemented

  int m_Imx;
  int m_Imy;
  int m_Imz;
protected:
  MABMISDisplacementFieldFilter();
  ~MABMISDisplacementFieldFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
};
} // namespace itk
} // namespace Statistics

#include "itkMABMISDisplacementFieldFilter.hxx"

#endif
