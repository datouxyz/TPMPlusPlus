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
#ifndef sitkBinaryErodeImageFilter_h
#define sitkBinaryErodeImageFilter_h

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

    /**\class BinaryErodeImageFilter
\brief Fast binary erosion of a single intensity value in the image.

BinaryErodeImageFilter is a binary erosion morphologic operation on the foreground of an image. Only the value designated by the intensity value "SetForegroundValue()" (alias as SetErodeValue() ) is considered as foreground, and other intensity values are considered background.

Gray scale images can be processed as binary images by selecting a "ForegroundValue" (alias "ErodeValue"). Pixel values matching the erode value are considered the "foreground" and all other pixels are "background". This is useful in processing segmented images where all pixels in segment #1 have value 1 and pixels in segment #2 have value 2, etc. A particular "segment number" can be processed. ForegroundValue defaults to the maximum possible value of the PixelType. The eroded pixels will receive the BackgroundValue (defaults to NumericTraits::NonpositiveMin() ).

The structuring element is assumed to be composed of binary values (zero or one). Only elements of the structuring element having values > 0 are candidates for affecting the center pixel. A reasonable choice of structuring element is itk::BinaryBallStructuringElement .

This implementation is based on the papers:

L.Vincent "Morphological transformations of binary images with
arbitrary structuring elements", and

N.Nikopoulos et al. "An efficient algorithm for 3d binary morphological transformations with 3d structuring elements for arbitrary size and shape". IEEE Transactions on Image Processing. Vol. 9. No. 3. 2000. pp. 283-286.

\see ImageToImageFilter BinaryDilateImageFilter BinaryMorphologyImageFilter
\sa itk::simple::BinaryErode for the procedural interface
\sa itk::BinaryErodeImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT BinaryErodeImageFilter : public ImageFilter {
    public:
      using Self = BinaryErodeImageFilter;

      /** Destructor */
      virtual ~BinaryErodeImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      BinaryErodeImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;
\

      /**
       * Set the radius of the kernel structuring element.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetKernelRadius ( std::vector<unsigned int> KernelRadius ) { this->m_KernelRadius = std::move(KernelRadius); return *this; }

      /** Set the values of the KernelRadius vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetKernelRadius( unsigned int value ) { this->m_KernelRadius = std::vector<unsigned int>(3, value); return *this; }

      /**
       * Get the radius of the kernel structuring element.
       */
      std::vector<unsigned int> GetKernelRadius() const { return this->m_KernelRadius; }\

      /**
       * Set the kernel or structuring element used for the morphology.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetKernelType ( KernelEnum KernelType ) { this->m_KernelType = KernelType; return *this; }

      /**
       * Get the kernel or structuring element used for the morphology.
       */
      KernelEnum GetKernelType() const { return this->m_KernelType; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( double BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       */
      double GetBackgroundValue() const { return this->m_BackgroundValue; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( double ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       */
      double GetForegroundValue() const { return this->m_ForegroundValue; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBoundaryToForeground ( bool BoundaryToForeground ) { this->m_BoundaryToForeground = BoundaryToForeground; return *this; }

      /** Set the value of BoundaryToForeground to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER BoundaryToForegroundOn() { return this->SetBoundaryToForeground(true); }
      SITK_RETURN_SELF_TYPE_HEADER BoundaryToForegroundOff() { return this->SetBoundaryToForeground(false); }

      /**
       */
      bool GetBoundaryToForeground() const { return this->m_BoundaryToForeground; }

      /** Name of this class */
      std::string GetName() const { return std::string ("BinaryErodeImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */

      Image Execute ( const Image& image1 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      /*  */
      std::vector<unsigned int>  m_KernelRadius{std::vector<uint32_t>(3, 1)};

      KernelEnum  m_KernelType{itk::simple::sitkBall};

      /*  */
      double  m_BackgroundValue{0.0};

      /*  */
      double  m_ForegroundValue{1.0};

      /*  */
      bool  m_BoundaryToForeground{true};


    };

    /**\
     * \brief Fast binary erosion of a single intensity value in the image.
     *
     * This function directly calls the execute method of BinaryErodeImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::BinaryErodeImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image BinaryErode ( const Image& image1, std::vector<unsigned int> kernelRadius = std::vector<uint32_t>(3, 1), KernelEnum kernelType = itk::simple::sitkBall, double backgroundValue = 0.0, double foregroundValue = 1.0, bool boundaryToForeground = true );

     /** @} */
  }
}
#endif
