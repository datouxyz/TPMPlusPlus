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
#ifndef sitkOrImageFilter_h
#define sitkOrImageFilter_h

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

    /**\class OrImageFilter
\brief Implements the OR bitwise operator pixel-wise between two images.

This class is templated over the types of the two input images and the type of the output image. Numeric conversions (castings) are done by the C++ defaults.

Since the bitwise OR operation is only defined in C++ for integer types, the images passed to this filter must comply with the requirement of using integer pixel type.

The total operation over one pixel will be

\code
output_pixel = static_cast<OutputPixelType>( input1_pixel | input2_pixel )

\endcode


Where "|" is the boolean OR operator in C++.
\sa itk::simple::Or for the procedural interface
\sa itk::OrImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT OrImageFilter : public ImageFilter {
    public:
      using Self = OrImageFilter;

      /** Destructor */
      virtual ~OrImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      OrImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;



      /** Name of this class */
      std::string GetName() const { return std::string ("OrImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input images */
#ifndef SWIG
      Image Execute ( Image&& image1, const Image& image2 );
#endif
      Image Execute ( const Image& image1, const Image& image2 );

      /** Execute the filter with an image and a constant */
      Image Execute ( const Image& image1, int constant );
#ifndef SWIG
      Image Execute ( Image&& image1, int constant );
#endif
      Image Execute ( int constant, const Image& image2 );


    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1, const Image& image2 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1, const Image& image2 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;
      using MemberFunction1Type = Image (Self::*)( int constant, const Image& image2 );
      template <class TImageType> Image ExecuteInternal ( int constant, const Image& image2 );
      friend struct detail::MemberFunctionAddressor<MemberFunction1Type>;
      std::unique_ptr<detail::MemberFunctionFactory<MemberFunction1Type> > m_MemberFactory1;

      using MemberFunction2Type = Image (Self::*)( const Image& image1, int constant );
      template <class TImageType> Image ExecuteInternal ( const Image& image1, int constant );
      friend struct detail::MemberFunctionAddressor<MemberFunction2Type>;
      std::unique_ptr<detail::MemberFunctionFactory<MemberFunction2Type> > m_MemberFactory2;



      bool m_InPlace{false};
    };

    /**\
     * \brief Implements the OR bitwise operator pixel-wise between two images.
     *
     * This function directly calls the execute method of OrImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::OrImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image Or ( Image&& image1, const Image& image2 );
#endif
     SITKBasicFilters_EXPORT Image Or ( const Image& image1, const Image& image2 );

     /** @} */
     SITKBasicFilters_EXPORT Image Or ( const Image& image1, int constant );
#ifndef  SWIG
     SITKBasicFilters_EXPORT Image Or ( Image&& image1, int constant );
#endif
     SITKBasicFilters_EXPORT Image Or ( int constant, const Image& image2 );
  }
}
#endif
