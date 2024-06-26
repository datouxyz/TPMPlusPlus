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
#ifndef sitkBinaryThresholdProjectionImageFilter_h
#define sitkBinaryThresholdProjectionImageFilter_h

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

    /**\class BinaryThresholdProjectionImageFilter
\brief BinaryThreshold projection.

This class was contributed to the Insight Journal by Gaetan Lehmann. the original paper can be found at https://hdl.handle.net/1926/164 

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


\see ProjectionImageFilter 


\see MedianProjectionImageFilter 


\see MeanProjectionImageFilter 


\see MeanProjectionImageFilter 


\see MaximumProjectionImageFilter 


\see MinimumProjectionImageFilter 


\see StandardDeviationProjectionImageFilter 


\see SumProjectionImageFilter
\sa itk::simple::BinaryThresholdProjection for the procedural interface
\sa itk::BinaryThresholdProjectionImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT BinaryThresholdProjectionImageFilter : public ImageFilter {
    public:
      using Self = BinaryThresholdProjectionImageFilter;

      /** Destructor */
      virtual ~BinaryThresholdProjectionImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      BinaryThresholdProjectionImageFilter();

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
       * Set/Get the input value consider as "threshold". Defaults to NumericTraits<InputPixelType>::max()
       */
      SITK_RETURN_SELF_TYPE_HEADER SetThresholdValue ( double ThresholdValue ) { this->m_ThresholdValue = ThresholdValue; return *this; }

      /**
       * Set/Get the input value consider as "threshold". Defaults to NumericTraits<InputPixelType>::max()
       */
      double GetThresholdValue() const { return this->m_ThresholdValue; }\

      /**
       * Set/Get the output value used as "foreground". Defaults to maximum value of PixelType.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( uint8_t ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       * Set/Get the output value used as "foreground". Defaults to maximum value of PixelType.
       */
      uint8_t GetForegroundValue() const { return this->m_ForegroundValue; }\

      /**
       * Set/Get the output value used as "background". Defaults to NumericTraits<PixelType>::NonpositiveMin() .
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( uint8_t BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Set/Get the output value used as "background". Defaults to NumericTraits<PixelType>::NonpositiveMin() .
       */
      uint8_t GetBackgroundValue() const { return this->m_BackgroundValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("BinaryThresholdProjectionImageFilter"); }

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

      double  m_ThresholdValue{0.0};

      uint8_t  m_ForegroundValue{1u};

      uint8_t  m_BackgroundValue{0u};


    };

    /**\
     * \brief BinaryThreshold projection.
     *
     * This function directly calls the execute method of BinaryThresholdProjectionImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::BinaryThresholdProjectionImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image BinaryThresholdProjection ( const Image& image1, unsigned int projectionDimension = 0u, double thresholdValue = 0.0, uint8_t foregroundValue = 1u, uint8_t backgroundValue = 0u );

     /** @} */
  }
}
#endif
