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
#ifndef sitkLabelMapOverlayImageFilter_h
#define sitkLabelMapOverlayImageFilter_h

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

    /**\class LabelMapOverlayImageFilter
\brief Apply a colormap to a label map and superimpose it on an image.

Apply a colormap to a label map and put it on top of the feature image. The feature image is typically the image from which the labeling was produced. Use the SetInput function to set the LabelMap , and the SetFeatureImage function to set the feature image.

The set of colors is a good selection of distinct colors. The opacity of the label map can be defined by the user. A background label produce a gray pixel with the same intensity than the input one.

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


This implementation was taken from the Insight Journal paper: https://hdl.handle.net/1926/584 or http://www.insight-journal.org/browse/publication/176 

\see LabelOverlayImageFilter , LabelOverlayFunctor 


\see LabelMapToRGBImageFilter , LabelMapToBinaryImageFilter , LabelMapToLabelImageFilter
\sa itk::simple::LabelMapOverlay for the procedural interface
\sa itk::LabelMapOverlayImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT LabelMapOverlayImageFilter : public ImageFilter {
    public:
      using Self = LabelMapOverlayImageFilter;

      /** Destructor */
      virtual ~LabelMapOverlayImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      LabelMapOverlayImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = LabelPixelIDTypeList;

\

      /**
       * Set/Get the opacity of the colored label image. The value must be between 0 and 1
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOpacity ( double Opacity ) { this->m_Opacity = Opacity; return *this; }

      /**
       * Set/Get the opacity of the colored label image. The value must be between 0 and 1
       */
      double GetOpacity() const { return this->m_Opacity; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetColormap ( std::vector<unsigned char> Colormap ) { this->m_Colormap = Colormap; return *this; }

      /**
       */
      std::vector<unsigned char> GetColormap() const { return this->m_Colormap; }

      /** Name of this class */
      std::string GetName() const { return std::string ("LabelMapOverlayImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */

      Image Execute ( const Image & labelMapImage, const Image & featureImage );


    private:
      /** Setup for member function dispatching */
      using MemberFunctionType = Image (Self::*)( const Image * labelMapImage, const Image * featureImage );

      friend struct detail::DualExecuteInternalAddressor<MemberFunctionType>;
      template <class TImageType1, class TImageType2> Image DualExecuteInternal ( const Image * labelMapImage, const Image * featureImage );


      std::unique_ptr<detail::DualMemberFunctionFactory<MemberFunctionType> > m_DualMemberFactory;



      /* Value assigned to pixels outside of the mask */
      double  m_Opacity{0.5};

      std::vector<unsigned char>  m_Colormap{std::vector<unsigned char>()};


    };

    /**\
     * \brief Apply a colormap to a label map and superimpose it on an image.
     *
     * This function directly calls the execute method of LabelMapOverlayImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::LabelMapOverlayImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image LabelMapOverlay ( const Image & labelMapImage, const Image & featureImage, double opacity = 0.5, std::vector<unsigned char> colormap = std::vector<unsigned char>() );

     /** @} */
  }
}
#endif