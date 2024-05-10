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
#ifndef sitkBinaryReconstructionByDilationImageFilter_h
#define sitkBinaryReconstructionByDilationImageFilter_h

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

    /**\class BinaryReconstructionByDilationImageFilter
\brief binary reconstruction by dilation of an image

Reconstruction by dilation operates on a "marker" image and a "mask" image, and is defined as the dilation of the marker image with respect to the mask image iterated until stability.

Geodesic morphology is described in Chapter 6.2 of Pierre Soille's book "Morphological Image Analysis: Principles and Applications", Second Edition, Springer, 2003.

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


This implementation was taken from the Insight Journal paper: https://hdl.handle.net/1926/584 or http://www.insight-journal.org/browse/publication/176 

\see MorphologyImageFilter , ReconstructionByDilationImageFilter , BinaryReconstructionByErosionImageFilter
\sa itk::simple::BinaryReconstructionByDilation for the procedural interface
\sa itk::BinaryReconstructionByDilationImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT BinaryReconstructionByDilationImageFilter : public ImageFilter {
    public:
      using Self = BinaryReconstructionByDilationImageFilter;

      /** Destructor */
      virtual ~BinaryReconstructionByDilationImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      BinaryReconstructionByDilationImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;
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
       * Set/Get the value used as "foreground" in the output image. Defaults to NumericTraits<PixelType>::max() .
       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( double ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       * Set/Get the value used as "foreground" in the output image. Defaults to NumericTraits<PixelType>::max() .
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
      bool GetFullyConnected() const { return this->m_FullyConnected; }

      /** Name of this class */
      std::string GetName() const { return std::string ("BinaryReconstructionByDilationImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */

      Image Execute ( const Image & markerImage, const Image & maskImage );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image * markerImage, const Image * maskImage );
      template <class TImageType> Image ExecuteInternal ( const Image * markerImage, const Image * maskImage );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      /*  */
      double  m_BackgroundValue{0.0};

      /*  */
      double  m_ForegroundValue{1.0};

      /*  */
      bool  m_FullyConnected{false};


    };

    /**\
     * \brief binary reconstruction by dilation of an image
     *
     * This function directly calls the execute method of BinaryReconstructionByDilationImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::BinaryReconstructionByDilationImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image BinaryReconstructionByDilation ( const Image & markerImage, const Image & maskImage, double backgroundValue = 0.0, double foregroundValue = 1.0, bool fullyConnected = false );

     /** @} */
  }
}
#endif
