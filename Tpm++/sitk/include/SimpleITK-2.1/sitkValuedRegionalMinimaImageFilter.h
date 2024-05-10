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
#ifndef sitkValuedRegionalMinimaImageFilter_h
#define sitkValuedRegionalMinimaImageFilter_h

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

    /**\class ValuedRegionalMinimaImageFilter
\brief Transforms the image so that any pixel that is not a regional minima is set to the maximum value for the pixel type. Pixels that are regional minima retain their value.

Regional minima are flat zones surrounded by pixels of higher value. A completely flat image will be marked as a regional minima by this filter.

This code was contributed in the Insight Journal paper: "Finding regional extrema - methods and performance" by Beare R., Lehmann G. https://hdl.handle.net/1926/153 http://www.insight-journal.org/browse/publication/65 

\author Richard Beare. Department of Medicine, Monash University, Melbourne, Australia.


\see ValuedRegionalMaximaImageFilter , ValuedRegionalExtremaImageFilter , 


\see HMinimaImageFilter
\sa itk::simple::ValuedRegionalMinima for the procedural interface
\sa itk::ValuedRegionalMinimaImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT ValuedRegionalMinimaImageFilter : public ImageFilter {
    public:
      using Self = ValuedRegionalMinimaImageFilter;

      /** Destructor */
      virtual ~ValuedRegionalMinimaImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      ValuedRegionalMinimaImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = ScalarPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetFullyConnected ( bool FullyConnected ) { this->m_FullyConnected = FullyConnected; return *this; }

      /** Set the value of FullyConnected to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER FullyConnectedOn() { return this->SetFullyConnected(true); }
      SITK_RETURN_SELF_TYPE_HEADER FullyConnectedOff() { return this->SetFullyConnected(false); }

      /**
       */
      bool GetFullyConnected() const { return this->m_FullyConnected; }
     /**
      *
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     bool GetFlat() const { return this->m_Flat; };


      /** Name of this class */
      std::string GetName() const { return std::string ("ValuedRegionalMinimaImageFilter"); }

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


      bool  m_FullyConnected{false};


      bool m_Flat{false};


    };

    /**\
     * \brief Transforms the image so that any pixel that is not a regional minima is set to the maximum value for the pixel type. Pixels that are regional minima retain their value.
     *
     * This function directly calls the execute method of ValuedRegionalMinimaImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::ValuedRegionalMinimaImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image ValuedRegionalMinima ( const Image& image1, bool fullyConnected = false );

     /** @} */
  }
}
#endif