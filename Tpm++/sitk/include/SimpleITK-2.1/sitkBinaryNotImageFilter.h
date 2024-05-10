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
#ifndef sitkBinaryNotImageFilter_h
#define sitkBinaryNotImageFilter_h

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

    /**\class BinaryNotImageFilter
\brief Implements the BinaryNot logical operator pixel-wise between two images.

This class is parameterized over the types of the two input images and the type of the output image. Numeric conversions (castings) are done by the C++ defaults.

The total operation over one pixel will be

output_pixel = static_cast<PixelType>( input1_pixel != input2_pixel )

Where "!=" is the equality operator in C++.

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


This implementation was taken from the Insight Journal paper: https://hdl.handle.net/1926/584 or http://www.insight-journal.org/browse/publication/176
\sa itk::simple::BinaryNot for the procedural interface
\sa itk::BinaryNotImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT BinaryNotImageFilter : public ImageFilter {
    public:
      using Self = BinaryNotImageFilter;

      /** Destructor */
      virtual ~BinaryNotImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      BinaryNotImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;
\

      /**
       * Set/Get the value in the image considered as "foreground". Defaults to maximum value of PixelType.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( double ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       * Set/Get the value in the image considered as "foreground". Defaults to maximum value of PixelType.
       */
      double GetForegroundValue() const { return this->m_ForegroundValue; }\

      /**
       * Set the value used as "background". Defaults to NumericTraits<PixelType>::NonpositiveMin() .
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( double BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Get the value used as "background". Defaults to NumericTraits<PixelType>::NonpositiveMin() .
       */
      double GetBackgroundValue() const { return this->m_BackgroundValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("BinaryNotImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */
#ifndef SWIG
      Image Execute ( Image&& image1 );
#endif
      Image Execute ( const Image& image1 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      double  m_ForegroundValue{1.0};

      double  m_BackgroundValue{0.0};


      bool m_InPlace{false};
    };

    /**\
     * \brief Implements the BinaryNot logical operator pixel-wise between two images.
     *
     * This function directly calls the execute method of BinaryNotImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::BinaryNotImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image BinaryNot ( Image&& image1, double foregroundValue = 1.0, double backgroundValue = 0.0 );
#endif
     SITKBasicFilters_EXPORT Image BinaryNot ( const Image& image1, double foregroundValue = 1.0, double backgroundValue = 0.0 );

     /** @} */
  }
}
#endif
