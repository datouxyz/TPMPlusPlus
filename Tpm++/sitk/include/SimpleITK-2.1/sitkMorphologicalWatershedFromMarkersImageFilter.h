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
#ifndef sitkMorphologicalWatershedFromMarkersImageFilter_h
#define sitkMorphologicalWatershedFromMarkersImageFilter_h

/*
 * WARNING: DO NOT EDIT THIS FILE!
 * THIS FILE IS AUTOMATICALLY GENERATED BY THE SIMPLEITK BUILD PROCESS.
 * Please look at sitkDualImageFilterTemplate.h.in to make changes.
 */

#include <memory>

#include "sitkImageFilter.h"
#include "sitkDualMemberFunctionFactory.h"
#include "sitkBasicFilters.h"

namespace itk {
  namespace simple {

    /**\class MorphologicalWatershedFromMarkersImageFilter
\brief Morphological watershed transform from markers.

The watershed transform is a tool for image segmentation that is fast and flexible and potentially fairly parameter free. It was originally derived from a geophysical model of rain falling on a terrain and a variety of more formal definitions have been devised to allow development of practical algorithms. If an image is considered as a terrain and divided into catchment basins then the hope is that each catchment basin would contain an object of interest.

The output is a label image. A label image, sometimes referred to as a categorical image, has unique values for each region. For example, if a watershed produces 2 regions, all pixels belonging to one region would have value A, and all belonging to the other might have value B. Unassigned pixels, such as watershed lines, might have the background value (0 by convention).

The simplest way of using the watershed is to preprocess the image we want to segment so that the boundaries of our objects are bright (e.g apply an edge detector) and compute the watershed transform of the edge image. Watershed lines will correspond to the boundaries and our problem will be solved. This is rarely useful in practice because there are always more regional minima than there are objects, either due to noise or natural variations in the object surfaces. Therefore, while many watershed lines do lie on significant boundaries, there are many that don't. Various methods can be used to reduce the number of minima in the image, like thresholding the smallest values, filtering the minima and/or smoothing the image.

This filter use another approach to avoid the problem of over segmentation: it let the user provide a marker image which mark the minima in the input image and give them a label. The minima are imposed in the input image by the markers. The labels of the output image are the label of the marker image.

The morphological watershed transform algorithm is described in Chapter 9.2 of Pierre Soille's book "Morphological Image Analysis:
 Principles and Applications", Second Edition, Springer, 2003.

This code was contributed in the Insight Journal paper: "The watershed transform in ITK - discussion and new developments" by Beare R., Lehmann G. https://hdl.handle.net/1926/202 http://www.insight-journal.org/browse/publication/92 

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France. 


\author Richard Beare. Department of Medicine, Monash University, Melbourne, Australia.


\see WatershedImageFilter , MorphologicalWatershedImageFilter
\sa itk::simple::MorphologicalWatershedFromMarkers for the procedural interface
\sa itk::MorphologicalWatershedFromMarkersImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT MorphologicalWatershedFromMarkersImageFilter : public ImageFilter {
    public:
      using Self = MorphologicalWatershedFromMarkersImageFilter;

      /** Destructor */
      virtual ~MorphologicalWatershedFromMarkersImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      MorphologicalWatershedFromMarkersImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = ScalarPixelIDTypeList;

\

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
      std::string GetName() const { return std::string ("MorphologicalWatershedFromMarkersImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */

      Image Execute ( const Image & image, const Image & markerImage );


    private:
      /** Setup for member function dispatching */
      using MemberFunctionType = Image (Self::*)( const Image * image, const Image * markerImage );

      friend struct detail::DualExecuteInternalAddressor<MemberFunctionType>;
      template <class TImageType1, class TImageType2> Image DualExecuteInternal ( const Image * image, const Image * markerImage );


      std::unique_ptr<detail::DualMemberFunctionFactory<MemberFunctionType> > m_DualMemberFactory;



      bool  m_MarkWatershedLine{true};

      bool  m_FullyConnected{false};


    };

    /**\
     * \brief Morphological watershed transform from markers.
     *
     * This function directly calls the execute method of MorphologicalWatershedFromMarkersImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::MorphologicalWatershedFromMarkersImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image MorphologicalWatershedFromMarkers ( const Image & image, const Image & markerImage, bool markWatershedLine = true, bool fullyConnected = false );

     /** @} */
  }
}
#endif