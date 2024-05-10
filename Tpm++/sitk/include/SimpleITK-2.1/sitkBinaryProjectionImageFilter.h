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
#ifndef sitkBinaryProjectionImageFilter_h
#define sitkBinaryProjectionImageFilter_h

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

    /**\class BinaryProjectionImageFilter
\brief Binary projection.

This class was contributed to the Insight Journal by Gaetan Lehmann. The original paper can be found at https://hdl.handle.net/1926/164 

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


\see ProjectionImageFilter 


\see MedianProjectionImageFilter 


\see MeanProjectionImageFilter 


\see MeanProjectionImageFilter 


\see MaximumProjectionImageFilter 


\see MinimumProjectionImageFilter 


\see StandardDeviationProjectionImageFilter 


\see SumProjectionImageFilter
\sa itk::simple::BinaryProjection for the procedural interface
\sa itk::BinaryProjectionImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT BinaryProjectionImageFilter : public ImageFilter {
    public:
      using Self = BinaryProjectionImageFilter;

      /** Destructor */
      virtual ~BinaryProjectionImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      BinaryProjectionImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetProjectionDimension ( unsigned int ProjectionDimension ) { this->m_ProjectionDimension = ProjectionDimension; return *this; }

      /**
       */
      unsigned int GetProjectionDimension() const { return this->m_ProjectionDimension; }\

      /**
       * Set the value in the image to consider as "foreground". Defaults to maximum value of PixelType. Subclasses may alias this to DilateValue or ErodeValue.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( double ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       * Get the value in the image considered as "foreground". Defaults to maximum value of PixelType.
       */
      double GetForegroundValue() const { return this->m_ForegroundValue; }\

      /**
       * Set the value used as "background". Any pixel value which is not DilateValue is considered background. BackgroundValue is used for defining boundary conditions. Defaults to NumericTraits<PixelType>::NonpositiveMin() .
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( double BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Get the value used as "background". Any pixel value which is not DilateValue is considered background. BackgroundValue is used for defining boundary conditions. Defaults to NumericTraits<PixelType>::NonpositiveMin() .
       */
      double GetBackgroundValue() const { return this->m_BackgroundValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("BinaryProjectionImageFilter"); }

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


      unsigned int  m_ProjectionDimension{0u};

      double  m_ForegroundValue{1.0};

      double  m_BackgroundValue{0.0};


    };

    /**\
     * \brief Binary projection.
     *
     * This function directly calls the execute method of BinaryProjectionImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::BinaryProjectionImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image BinaryProjection ( const Image& image1, unsigned int projectionDimension = 0u, double foregroundValue = 1.0, double backgroundValue = 0.0 );

     /** @} */
  }
}
#endif