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
#ifndef sitkTernaryMagnitudeSquaredImageFilter_h
#define sitkTernaryMagnitudeSquaredImageFilter_h

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

    /**\class TernaryMagnitudeSquaredImageFilter
\brief Compute the pixel-wise squared magnitude of three images.

This class is templated over the types of the three input images and the type of the output image. Numeric conversions (castings) are done by the C++ defaults.
\sa itk::simple::TernaryMagnitudeSquared for the procedural interface
\sa itk::TernaryMagnitudeSquaredImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT TernaryMagnitudeSquaredImageFilter : public ImageFilter {
    public:
      using Self = TernaryMagnitudeSquaredImageFilter;

      /** Destructor */
      virtual ~TernaryMagnitudeSquaredImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      TernaryMagnitudeSquaredImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;


      /** Name of this class */
      std::string GetName() const { return std::string ("TernaryMagnitudeSquaredImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */
#ifndef SWIG
      Image Execute ( Image&& image1, const Image& image2, const Image& image3 );
#endif
      Image Execute ( const Image& image1, const Image& image2, const Image& image3 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1, const Image& image2, const Image& image3 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1, const Image& image2, const Image& image3 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;



      bool m_InPlace{false};
    };

    /**\
     * \brief Compute the pixel-wise squared magnitude of three images.
     *
     * This function directly calls the execute method of TernaryMagnitudeSquaredImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::TernaryMagnitudeSquaredImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image TernaryMagnitudeSquared ( Image&& image1, const Image& image2, const Image& image3 );
#endif
     SITKBasicFilters_EXPORT Image TernaryMagnitudeSquared ( const Image& image1, const Image& image2, const Image& image3 );

     /** @} */
  }
}
#endif
