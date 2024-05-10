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
#ifndef sitkDivideFloorImageFilter_h
#define sitkDivideFloorImageFilter_h

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

    /**\class DivideFloorImageFilter
\brief Implements pixel-wise generic operation of two images, or of an image and a constant.

This class is parameterized over the types of the two input images and the type of the output image. It is also parameterized by the operation to be applied. A Functor style is used.

The constant must be of the same type than the pixel type of the corresponding image. It is wrapped in a SimpleDataObjectDecorator so it can be updated through the pipeline. The SetConstant() and GetConstant() methods are provided as shortcuts to set or get the constant value without manipulating the decorator.

\see BinaryGeneratorImagFilter 


\see UnaryFunctorImageFilter TernaryFunctorImageFilter
\sa itk::simple::DivideFloor for the procedural interface
\sa itk::BinaryFunctorImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT DivideFloorImageFilter : public ImageFilter {
    public:
      using Self = DivideFloorImageFilter;

      /** Destructor */
      virtual ~DivideFloorImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      DivideFloorImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;



      /** Name of this class */
      std::string GetName() const { return std::string ("DivideFloorImageFilter"); }

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
     * \brief Implements pixel-wise generic operation of two images, or of an image and a constant.
     *
     * This function directly calls the execute method of DivideFloorImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::DivideFloorImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image DivideFloor ( Image&& image1, const Image& image2 );
#endif
     SITKBasicFilters_EXPORT Image DivideFloor ( const Image& image1, const Image& image2 );

     /** @} */
     SITKBasicFilters_EXPORT Image DivideFloor ( const Image& image1, double constant );
#ifndef  SWIG
     SITKBasicFilters_EXPORT Image DivideFloor ( Image&& image1, double constant );
#endif
     SITKBasicFilters_EXPORT Image DivideFloor ( double constant, const Image& image2 );
  }
}
#endif
