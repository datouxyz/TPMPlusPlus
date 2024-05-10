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
#ifndef sitkLessEqualImageFilter_h
#define sitkLessEqualImageFilter_h

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

    /**\class LessEqualImageFilter
\brief Implements pixel-wise generic operation of two images, or of an image and a constant.

This class is parameterized over the types of the two input images and the type of the output image. It is also parameterized by the operation to be applied. A Functor style is used.

The constant must be of the same type than the pixel type of the corresponding image. It is wrapped in a SimpleDataObjectDecorator so it can be updated through the pipeline. The SetConstant() and GetConstant() methods are provided as shortcuts to set or get the constant value without manipulating the decorator.

\see BinaryGeneratorImagFilter 


\see UnaryFunctorImageFilter TernaryFunctorImageFilter
\sa itk::simple::LessEqual for the procedural interface
\sa itk::BinaryFunctorImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT LessEqualImageFilter : public ImageFilter {
    public:
      using Self = LessEqualImageFilter;

      /** Destructor */
      virtual ~LessEqualImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      LessEqualImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;

\

      /**
       * Set/Get the value used to mark the false pixels of the operator.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( uint8_t BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Set/Get the value used to mark the false pixels of the operator.

       */
      uint8_t GetBackgroundValue() const { return this->m_BackgroundValue; }\

      /**
       * Set/Get the value used to mark the true pixels of the operator.

       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( uint8_t ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       * Set/Get the value used to mark the true pixels of the operator.

       */
      uint8_t GetForegroundValue() const { return this->m_ForegroundValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("LessEqualImageFilter"); }

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
      /** Execute the filter on an image and a constant with the given parameters */
      Image Execute ( const Image& image1, double constant, uint8_t backgroundValue, uint8_t foregroundValue );
      Image Execute ( double constant, const Image& image2, uint8_t backgroundValue, uint8_t foregroundValue );

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


      uint8_t  m_BackgroundValue{0u};

      uint8_t  m_ForegroundValue{1u};


      bool m_InPlace{false};
    };

    /**\
     * \brief Implements pixel-wise generic operation of two images, or of an image and a constant.
     *
     * This function directly calls the execute method of LessEqualImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::LessEqualImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image LessEqual ( Image&& image1, const Image& image2, uint8_t backgroundValue = 0u, uint8_t foregroundValue = 1u );
#endif
     SITKBasicFilters_EXPORT Image LessEqual ( const Image& image1, const Image& image2, uint8_t backgroundValue = 0u, uint8_t foregroundValue = 1u );

     /** @} */
     SITKBasicFilters_EXPORT Image LessEqual ( const Image& image1, double constant, uint8_t backgroundValue = 0u, uint8_t foregroundValue = 1u );
#ifndef  SWIG
     SITKBasicFilters_EXPORT Image LessEqual ( Image&& image1, double constant, uint8_t backgroundValue = 0u, uint8_t foregroundValue = 1u );
#endif
     SITKBasicFilters_EXPORT Image LessEqual ( double constant, const Image& image2, uint8_t backgroundValue = 0u, uint8_t foregroundValue = 1u );
  }
}
#endif
