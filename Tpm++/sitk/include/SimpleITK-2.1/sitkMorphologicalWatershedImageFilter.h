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
#ifndef sitkMorphologicalWatershedImageFilter_h
#define sitkMorphologicalWatershedImageFilter_h

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

    /**\class MorphologicalWatershedImageFilter
\brief Watershed segmentation implementation with morphological operators.

Watershed pixel are labeled 0. TOutputImage should be an integer type. Labels of output image are in no particular order. You can reorder the labels such that object labels are consecutive and sorted based on object size by passing the output of this filter to a RelabelComponentImageFilter .

The morphological watershed transform algorithm is described in Chapter 9.2 of Pierre Soille's book "Morphological Image Analysis:
Principles and Applications", Second Edition, Springer, 2003.

This code was contributed in the Insight Journal paper: "The watershed transform in ITK - discussion and new developments" by Beare R., Lehmann G. https://hdl.handle.net/1926/202 http://www.insight-journal.org/browse/publication/92 

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


\see WatershedImageFilter , MorphologicalWatershedFromMarkersImageFilter
\sa itk::simple::MorphologicalWatershed for the procedural interface
\sa itk::MorphologicalWatershedImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT MorphologicalWatershedImageFilter : public ImageFilter {
    public:
      using Self = MorphologicalWatershedImageFilter;

      /** Destructor */
      virtual ~MorphologicalWatershedImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      MorphologicalWatershedImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = ScalarPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetLevel ( double Level ) { this->m_Level = Level; return *this; }

      /**
       */
      double GetLevel() const { return this->m_Level; }\

      /**
       * Set/Get whether the watershed pixel must be marked or not. Default is true. Set it to false do not only avoid writing watershed pixels, it also decrease algorithm complexity.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMarkWatershedLine ( bool MarkWatershedLine ) { this->m_MarkWatershedLine = MarkWatershedLine; return *this; }

      /** Set the value of MarkWatershedLine to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER MarkWatershedLineOn() { return this->SetMarkWatershedLine(true); }
      SITK_RETURN_SELF_TYPE_HEADER MarkWatershedLineOff() { return this->SetMarkWatershedLine(false); }

      /**
       * Set/Get whether the watershed pixel must be marked or not. Default is true. Set it to false do not only avoid writing watershed pixels, it also decrease algorithm complexity.
       */
      bool GetMarkWatershedLine() const { return this->m_MarkWatershedLine; }\

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
      bool GetFullyConnected() const { return this->m_FullyConnected; }

      /** Name of this class */
      std::string GetName() const { return std::string ("MorphologicalWatershedImageFilter"); }

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


      double  m_Level{0.0};

      bool  m_MarkWatershedLine{true};

      /*  */
      bool  m_FullyConnected{false};


    };

    /**\
     * \brief Watershed segmentation implementation with morphological operators.
     *
     * This function directly calls the execute method of MorphologicalWatershedImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::MorphologicalWatershedImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image MorphologicalWatershed ( const Image& image1, double level = 0.0, bool markWatershedLine = true, bool fullyConnected = false );

     /** @} */
  }
}
#endif
