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
#ifndef sitkRegionalMaximaImageFilter_h
#define sitkRegionalMaximaImageFilter_h

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

    /**\class RegionalMaximaImageFilter
\brief Produce a binary image where foreground is the regional maxima of the input image.

Regional maxima are flat zones surrounded by pixels of lower value.

If the input image is constant, the entire image can be considered as a maxima or not. The desired behavior can be selected with the SetFlatIsMaxima() method.

\author Gaetan Lehmann


This class was contributed to the Insight Journal by author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France. The paper can be found at https://hdl.handle.net/1926/153 

\see ValuedRegionalMaximaImageFilter 


\see HConvexImageFilter 


\see RegionalMinimaImageFilter
\sa itk::simple::RegionalMaxima for the procedural interface
\sa itk::RegionalMaximaImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT RegionalMaximaImageFilter : public ImageFilter {
    public:
      using Self = RegionalMaximaImageFilter;

      /** Destructor */
      virtual ~RegionalMaximaImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      RegionalMaximaImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = ScalarPixelIDTypeList;
\

      /**
       * Set/Get the value used as "background" in the output image. Defaults to NumericTraits<PixelType>::NonpositiveMin() .
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( double BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Set/Get the value used as "background" in the output image. Defaults to NumericTraits<PixelType>::NonpositiveMin() .
       */
      double GetBackgroundValue() const { return this->m_BackgroundValue; }\

      /**
       * Set/Get the value in the output image to consider as "foreground". Defaults to maximum value of PixelType.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( double ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       * Set/Get the value in the output image to consider as "foreground". Defaults to maximum value of PixelType.
       */
      double GetForegroundValue() const { return this->m_ForegroundValue; }\

      /**
       * Set/Get whether the connected components are defined strictly by face connectivity or by face+edge+vertex connectivity. Default is FullyConnectedOff. For objects that are 1 pixel wide, use FullyConnectedOn.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetFullyConnected ( bool FullyConnected ) { this->m_FullyConnected = FullyConnected; return *this; }

      /** Set the value of FullyConnected to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER FullyConnectedOn() { return this->SetFullyConnected(true); }
      SITK_RETURN_SELF_TYPE_HEADER FullyConnectedOff() { return this->SetFullyConnected(false); }

      /**
       * Set/Get whether the connected components are defined strictly by face connectivity or by face+edge+vertex connectivity. Default is FullyConnectedOff. For objects that are 1 pixel wide, use FullyConnectedOn.
       */
      bool GetFullyConnected() const { return this->m_FullyConnected; }\

      /**
       * Set/Get whether a flat image must be considered as a maxima or not. Defaults to true.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetFlatIsMaxima ( bool FlatIsMaxima ) { this->m_FlatIsMaxima = FlatIsMaxima; return *this; }

      /** Set the value of FlatIsMaxima to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER FlatIsMaximaOn() { return this->SetFlatIsMaxima(true); }
      SITK_RETURN_SELF_TYPE_HEADER FlatIsMaximaOff() { return this->SetFlatIsMaxima(false); }

      /**
       * Set/Get wether a flat image must be considered as a maxima or not. Defaults to true.
       */
      bool GetFlatIsMaxima() const { return this->m_FlatIsMaxima; }

      /** Name of this class */
      std::string GetName() const { return std::string ("RegionalMaximaImageFilter"); }

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
      double  m_BackgroundValue{0.0};

      /*  */
      double  m_ForegroundValue{1.0};

      /*  */
      bool  m_FullyConnected{false};

      bool  m_FlatIsMaxima{true};


    };

    /**\
     * \brief Produce a binary image where foreground is the regional maxima of the input image.
     *
     * This function directly calls the execute method of RegionalMaximaImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::RegionalMaximaImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image RegionalMaxima ( const Image& image1, double backgroundValue = 0.0, double foregroundValue = 1.0, bool fullyConnected = false, bool flatIsMaxima = true );

     /** @} */
  }
}
#endif
