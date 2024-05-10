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
#ifndef sitkBinaryMorphologicalOpeningImageFilter_h
#define sitkBinaryMorphologicalOpeningImageFilter_h

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

    /**\class BinaryMorphologicalOpeningImageFilter
\brief binary morphological opening of an image.

This filter removes small (i.e., smaller than the structuring element) structures in the interior or at the boundaries of the image. The morphological opening of an image "f" is defined as: Opening(f) = Dilatation(Erosion(f)).

The structuring element is assumed to be composed of binary values (zero or one). Only elements of the structuring element having values > 0 are candidates for affecting the center pixel.

This code was contributed in the Insight Journal paper: "Binary morphological closing and opening image filters" by Lehmann G. https://hdl.handle.net/1926/141 http://www.insight-journal.org/browse/publication/58 

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


\see MorphologyImageFilter , GrayscaleDilateImageFilter , GrayscaleErodeImageFilter
\sa itk::simple::BinaryMorphologicalOpening for the procedural interface
\sa itk::BinaryMorphologicalOpeningImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT BinaryMorphologicalOpeningImageFilter : public ImageFilter {
    public:
      using Self = BinaryMorphologicalOpeningImageFilter;

      /** Destructor */
      virtual ~BinaryMorphologicalOpeningImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      BinaryMorphologicalOpeningImageFilter();

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
       * Set the value in eroded part of the image. Defaults to zero
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( double BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Set the value in eroded part of the image. Defaults to zero
       */
      double GetBackgroundValue() const { return this->m_BackgroundValue; }\

      /**
       * Set the value in the image to consider as "foreground". Defaults to maximum value of PixelType.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( double ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       * Get the value in the image considered as "foreground". Defaults to maximum value of PixelType.
       */
      double GetForegroundValue() const { return this->m_ForegroundValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("BinaryMorphologicalOpeningImageFilter"); }

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

      double  m_BackgroundValue{0.0};

      double  m_ForegroundValue{1.0};


    };

    /**\
     * \brief binary morphological opening of an image.
     *
     * This function directly calls the execute method of BinaryMorphologicalOpeningImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::BinaryMorphologicalOpeningImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image BinaryMorphologicalOpening ( const Image& image1, std::vector<unsigned int> kernelRadius = std::vector<uint32_t>(3, 1), KernelEnum kernelType = itk::simple::sitkBall, double backgroundValue = 0.0, double foregroundValue = 1.0 );

     /** @} */
  }
}
#endif
