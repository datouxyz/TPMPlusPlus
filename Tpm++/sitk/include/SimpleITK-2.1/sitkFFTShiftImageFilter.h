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
#ifndef sitkFFTShiftImageFilter_h
#define sitkFFTShiftImageFilter_h

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

    /**\class FFTShiftImageFilter
\brief Shift the zero-frequency components of a Fourier transform to the center of the image.

The Fourier transform produces an image where the zero frequency components are in the corner of the image, making it difficult to understand. This filter shifts the component to the center of the image.

\note For images with an odd-sized dimension, applying this filter twice will not produce the same image as the original one without using SetInverse(true) on one (and only one) of the two filters.


https://hdl.handle.net/1926/321 

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


\see ForwardFFTImageFilter , InverseFFTImageFilter
\sa itk::simple::FFTShift for the procedural interface
\sa itk::FFTShiftImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT FFTShiftImageFilter : public ImageFilter {
    public:
      using Self = FFTShiftImageFilter;

      /** Destructor */
      virtual ~FFTShiftImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      FFTShiftImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = NonLabelPixelIDTypeList;
\

      /**
       * Set/Get whether the filter must invert the transform or not. This option has no effect if none of the size of the input image is even, but is required to restore the original image if at least one of the dimensions has an odd size.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetInverse ( bool Inverse ) { this->m_Inverse = Inverse; return *this; }

      /** Set the value of Inverse to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER InverseOn() { return this->SetInverse(true); }
      SITK_RETURN_SELF_TYPE_HEADER InverseOff() { return this->SetInverse(false); }

      /**
       * Set/Get whether the filter must invert the transform or not. This option has no effect if none of the size of the input image is even, but is required to restore the original image if at least one of the dimensions has an odd size.
       */
      bool GetInverse() const { return this->m_Inverse; }

      /** Name of this class */
      std::string GetName() const { return std::string ("FFTShiftImageFilter"); }

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
      bool  m_Inverse{false};


    };

    /**\
     * \brief Shift the zero-frequency components of a Fourier transform to the center of the image.
     *
     * This function directly calls the execute method of FFTShiftImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::FFTShiftImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image FFTShift ( const Image& image1, bool inverse = false );

     /** @} */
  }
}
#endif
