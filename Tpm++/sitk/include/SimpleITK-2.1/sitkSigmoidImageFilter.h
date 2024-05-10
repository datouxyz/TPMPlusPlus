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
#ifndef sitkSigmoidImageFilter_h
#define sitkSigmoidImageFilter_h

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

    /**\class SigmoidImageFilter
\brief Computes the sigmoid function pixel-wise.

A linear transformation is applied first on the argument of the sigmoid function. The resulting total transform is given by

 \f[ f(x) = (Max-Min) \cdot \frac{1}{\left(1+e^{- \frac{ x - \beta }{\alpha}}\right)} + Min \f] 

Every output pixel is equal to f(x). Where x is the intensity of the homologous input pixel, and alpha and beta are user-provided constants.
\sa itk::simple::Sigmoid for the procedural interface
\sa itk::SigmoidImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT SigmoidImageFilter : public ImageFilter {
    public:
      using Self = SigmoidImageFilter;

      /** Destructor */
      virtual ~SigmoidImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      SigmoidImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetAlpha ( double Alpha ) { this->m_Alpha = Alpha; return *this; }

      /**
       */
      double GetAlpha() const { return this->m_Alpha; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBeta ( double Beta ) { this->m_Beta = Beta; return *this; }

      /**
       */
      double GetBeta() const { return this->m_Beta; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOutputMaximum ( double OutputMaximum ) { this->m_OutputMaximum = OutputMaximum; return *this; }

      /**
       */
      double GetOutputMaximum() const { return this->m_OutputMaximum; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOutputMinimum ( double OutputMinimum ) { this->m_OutputMinimum = OutputMinimum; return *this; }

      /**
       */
      double GetOutputMinimum() const { return this->m_OutputMinimum; }

      /** Name of this class */
      std::string GetName() const { return std::string ("SigmoidImageFilter"); }

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


      /* Alpha */
      double  m_Alpha{1};

      /* Alpha */
      double  m_Beta{0};

      double  m_OutputMaximum{255};

      double  m_OutputMinimum{0};


      bool m_InPlace{false};
    };

    /**\
     * \brief Computes the sigmoid function pixel-wise.
     *
     * This function directly calls the execute method of SigmoidImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::SigmoidImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image Sigmoid ( Image&& image1, double alpha = 1, double beta = 0, double outputMaximum = 255, double outputMinimum = 0 );
#endif
     SITKBasicFilters_EXPORT Image Sigmoid ( const Image& image1, double alpha = 1, double beta = 0, double outputMaximum = 255, double outputMinimum = 0 );

     /** @} */
  }
}
#endif