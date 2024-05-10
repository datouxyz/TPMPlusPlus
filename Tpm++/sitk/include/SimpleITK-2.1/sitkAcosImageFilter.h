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
#ifndef sitkAcosImageFilter_h
#define sitkAcosImageFilter_h

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

    /**\class AcosImageFilter
\brief Computes the inverse cosine of each pixel.

This filter is templated over the pixel type of the input image and the pixel type of the output image.

The filter walks over all the pixels in the input image, and for each pixel does do the following:



\li cast the pixel value to double , 


\li apply the std::acos() function to the double value 


\li cast the double value resulting from std::acos() to the pixel type of the output image 


\li store the casted value into the output image.



The filter expects both images to have the same dimension (e.g. both 2D, or both 3D, or both ND).
\sa itk::simple::Acos for the procedural interface
\sa itk::AcosImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT AcosImageFilter : public ImageFilter {
    public:
      using Self = AcosImageFilter;

      /** Destructor */
      virtual ~AcosImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      AcosImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;


      /** Name of this class */
      std::string GetName() const { return std::string ("AcosImageFilter"); }

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
      /** Dispatched methods which calls ExecuteInteral on each component */
      template <class TImageType> Image ExecuteInternalVectorImage ( const Image& image );

      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;
      friend struct detail::ExecuteInternalVectorImageAddressor<MemberFunctionType>;
      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;



      bool m_InPlace{false};
    };

    /**\
     * \brief Computes the inverse cosine of each pixel.
     *
     * This function directly calls the execute method of AcosImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::AcosImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image Acos ( Image&& image1 );
#endif
     SITKBasicFilters_EXPORT Image Acos ( const Image& image1 );

     /** @} */
  }
}
#endif