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
#ifndef sitkPowImageFilter_h
#define sitkPowImageFilter_h

/*
 * WARNING: DO NOT EDIT THIS FILE!
 * THIS FILE IS AUTOMATICALLY GENERATED BY THE SIMPLEITK BUILD PROCESS.
 * Please look at sitkBinaryFunctorFilterTemplate.h.in to make changes.
 */

#include <memory>

#include "sitkBasicFilters.h"
#include "sitkImageFilter.h"

namespace itk {
  namespace simple {

    /**\class PowImageFilter
\brief Computes the powers of 2 images.

This class is templated over the types of the two input images and the type of the output image. Numeric conversions (castings) are done by the C++ defaults.

The output of the pow function will be cast to the pixel type of the output image.

The total operation over one pixel will be \code
output_pixel = static_cast< TOutput >( std::pow(static_cast<RealType>(A),static_cast<RealType>(B)) );

\endcode


The pow function can be applied to two images with the following: \code
SetInput1( image1 );

SetInput2( image2 );

\endcode


Additionally, this filter can be used to raise every pixel of an image to a power of a constant by using \code
SetInput1( image1 );

SetConstant2( constant );

\endcode
\sa itk::simple::Pow for the procedural interface
\sa itk::PowImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT PowImageFilter : public ImageFilter {
    public:
      using Self = PowImageFilter;

      /** Destructor */
      virtual ~PowImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      PowImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = typelist::Append<BasicPixelIDTypeList, ComplexPixelIDTypeList>::Type;



      /** Name of this class */
      std::string GetName() const { return std::string ("PowImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input images */
#ifndef SWIG
      Image Execute ( Image&& image1, const Image& image2 );
#endif
      Image Execute ( const Image& image1, const Image& image2 );

      /** Execute the filter with an image and a constant */
      Image Execute ( const Image& image1, double constant );
#ifndef SWIG
      Image Execute ( Image&& image1, double constant );
#endif
      Image Execute ( double constant, const Image& image2 );


    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1, const Image& image2 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1, const Image& image2 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;
      using MemberFunction1Type = Image (Self::*)( double constant, const Image& image2 );
      template <class TImageType> Image ExecuteInternal ( double constant, const Image& image2 );
      friend struct detail::MemberFunctionAddressor<MemberFunction1Type>;
      std::unique_ptr<detail::MemberFunctionFactory<MemberFunction1Type> > m_MemberFactory1;

      using MemberFunction2Type = Image (Self::*)( const Image& image1, double constant );
      template <class TImageType> Image ExecuteInternal ( const Image& image1, double constant );
      friend struct detail::MemberFunctionAddressor<MemberFunction2Type>;
      std::unique_ptr<detail::MemberFunctionFactory<MemberFunction2Type> > m_MemberFactory2;



      bool m_InPlace{false};
    };

    /**\
     * \brief Computes the powers of 2 images.
     *
     * This function directly calls the execute method of PowImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::PowImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image Pow ( Image&& image1, const Image& image2 );
#endif
     SITKBasicFilters_EXPORT Image Pow ( const Image& image1, const Image& image2 );

     /** @} */
     SITKBasicFilters_EXPORT Image Pow ( const Image& image1, double constant );
#ifndef  SWIG
     SITKBasicFilters_EXPORT Image Pow ( Image&& image1, double constant );
#endif
     SITKBasicFilters_EXPORT Image Pow ( double constant, const Image& image2 );
  }
}
#endif