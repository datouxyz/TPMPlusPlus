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
#ifndef sitkZeroCrossingImageFilter_h
#define sitkZeroCrossingImageFilter_h

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

    /**\class ZeroCrossingImageFilter
\brief This filter finds the closest pixel to the zero-crossings (sign changes) in a signed itk::Image .

Pixels closest to zero-crossings are labeled with a foreground value. All other pixels are marked with a background value. The algorithm works by detecting differences in sign among neighbors using city-block style connectivity (4-neighbors in 2d, 6-neighbors in 3d, etc.).

\par Inputs and Outputs
The input to this filter is an itk::Image of arbitrary dimension. The algorithm assumes a signed data type (zero-crossings are not defined for unsigned data types), and requires that operator>, operator<, operator==, and operator!= are defined.


\par 
The output of the filter is a binary, labeled image of user-specified type. By default, zero-crossing pixels are labeled with a default "foreground" value of itk::NumericTraits<OutputDataType>::OneValue() , where OutputDataType is the data type of the output image. All other pixels are labeled with a default "background" value of itk::NumericTraits<OutputDataType>::ZeroValue() .


\par Parameters
There are two parameters for this filter. ForegroundValue is the value that marks zero-crossing pixels. The BackgroundValue is the value given to all other pixels.


\see Image 


\see Neighborhood 


\see NeighborhoodOperator 


\see NeighborhoodIterator
\sa itk::simple::ZeroCrossing for the procedural interface
\sa itk::ZeroCrossingImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT ZeroCrossingImageFilter : public ImageFilter {
    public:
      using Self = ZeroCrossingImageFilter;

      /** Destructor */
      virtual ~ZeroCrossingImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      ZeroCrossingImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = SignedPixelIDTypeList;
\

      /**
       * Set/Get the label value for zero-crossing pixels.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( uint8_t ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       * Set/Get the label value for zero-crossing pixels.
       */
      uint8_t GetForegroundValue() const { return this->m_ForegroundValue; }\

      /**
       * Set/Get the label value for non-zero-crossing pixels.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( uint8_t BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Set/Get the label value for non-zero-crossing pixels.
       */
      uint8_t GetBackgroundValue() const { return this->m_BackgroundValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("ZeroCrossingImageFilter"); }

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


      uint8_t  m_ForegroundValue{1u};

      uint8_t  m_BackgroundValue{0u};


    };

    /**\
     * \brief This filter finds the closest pixel to the zero-crossings (sign changes) in a signed itk::Image .
     *
     * This function directly calls the execute method of ZeroCrossingImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::ZeroCrossingImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image ZeroCrossing ( const Image& image1, uint8_t foregroundValue = 1u, uint8_t backgroundValue = 0u );

     /** @} */
  }
}
#endif
