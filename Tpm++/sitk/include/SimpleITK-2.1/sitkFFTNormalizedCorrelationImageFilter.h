/*=========================================================================
*
*  Copyright NumFOCUS
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*         http://www.apache.org/licenses/LICENSE-2.0.txt
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*=========================================================================*/
#ifndef sitkFFTNormalizedCorrelationImageFilter_h
#define sitkFFTNormalizedCorrelationImageFilter_h

/*
 * WARNING: DO NOT EDIT THIS FILE!
 * THIS FILE IS AUTOMATICALLY GENERATED BY THE SIMPLEITK BUILD PROCESS.
 * Please look at sitkImageFilterTemplate.h.in to make changes.
 */

#include <memory>

#include "sitkBasicFilters.h"
#include "sitkImageFilter.h"

namespace itk {
  namespace simple {

    /**\class FFTNormalizedCorrelationImageFilter
\brief Calculate normalized cross correlation using FFTs.

This filter calculates the normalized cross correlation (NCC) of two images using FFTs instead of spatial correlation. It is much faster than spatial correlation for reasonably large structuring elements. This filter is a subclass of the more general MaskedFFTNormalizedCorrelationImageFilter and operates by essentially setting the masks in that algorithm to images of ones. As described in detail in the references below, there is no computational overhead to utilizing the more general masked algorithm because the FFTs of the images of ones are still necessary for the computations.

Inputs: Two images are required as inputs, fixedImage and movingImage. In the context of correlation, inputs are often defined as: "image" and "template". In this filter, the fixedImage plays the role of the image, and the movingImage plays the role of the template. However, this filter is capable of correlating any two images and is not restricted to small movingImages (templates).

Optional parameters: The RequiredNumberOfOverlappingPixels enables the user to specify how many voxels of the two images must overlap; any location in the correlation map that results from fewer than this number of voxels will be set to zero. Larger values zero-out pixels on a larger border around the correlation image. Thus, larger values remove less stable computations but also limit the capture range. If RequiredNumberOfOverlappingPixels is set to 0, the default, no zeroing will take place.

Image size: fixedImage and movingImage need not be the same size. Furthermore, whereas some algorithms require that the "template" be smaller than the "image" because of errors in the regions where the two are not fully overlapping, this filter has no such restriction.

Image spacing: Since the computations are done in the pixel domain, all input images must have the same spacing.

Outputs; The output is an image of RealPixelType that is the NCC of the two images and its values range from -1.0 to 1.0. The size of this NCC image is, by definition, size(fixedImage) + size(movingImage) - 1.

Example filter usage: \code
using FilterType = itk::FFTNormalizedCorrelationImageFilter< ShortImageType, DoubleImageType >;

FilterType::Pointer filter = FilterType::New();

filter->SetFixedImage( fixedImage );

filter->SetMovingImage( movingImage );

filter->SetRequiredNumberOfOverlappingPixels(20);

filter->Update();

\endcode


\warning The pixel type of the output image must be of real type (float or double). ConceptChecking is used to enforce the output pixel type. You will get a compilation error if the pixel type of the output image is not float or double.


References: 1) D. Padfield. "Masked object registration in the Fourier domain." Transactions on Image Processing. 2) D. Padfield. "Masked FFT registration". In Proc. Computer Vision and Pattern Recognition, 2010.

\author : Dirk Padfield, GE Global Research, padfield@research.ge.com
\sa itk::simple::FFTNormalizedCorrelation for the procedural interface
\sa itk::FFTNormalizedCorrelationImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT FFTNormalizedCorrelationImageFilter : public ImageFilter {
    public:
      using Self = FFTNormalizedCorrelationImageFilter;

      /** Destructor */
      virtual ~FFTNormalizedCorrelationImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      FFTNormalizedCorrelationImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetRequiredNumberOfOverlappingPixels ( uint64_t RequiredNumberOfOverlappingPixels ) { this->m_RequiredNumberOfOverlappingPixels = RequiredNumberOfOverlappingPixels; return *this; }

      /**
       */
      uint64_t GetRequiredNumberOfOverlappingPixels() const { return this->m_RequiredNumberOfOverlappingPixels; }

      /** Name of this class */
      std::string GetName() const { return std::string ("FFTNormalizedCorrelationImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input images */

      Image Execute ( const Image& image1, const Image& image2 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1, const Image& image2 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1, const Image& image2 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      uint64_t  m_RequiredNumberOfOverlappingPixels{0u};


    };

    /**\
     * \brief Calculate normalized cross correlation using FFTs.
     *
     * This function directly calls the execute method of FFTNormalizedCorrelationImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::FFTNormalizedCorrelationImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image FFTNormalizedCorrelation ( const Image& image1, const Image& image2, uint64_t requiredNumberOfOverlappingPixels = 0u );

     /** @} */
  }
}
#endif
