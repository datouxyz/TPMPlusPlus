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
#ifndef sitkLabelMapMaskImageFilter_h
#define sitkLabelMapMaskImageFilter_h

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

    /**\class LabelMapMaskImageFilter
\brief Mask and image with a LabelMap .

LabelMapMaskImageFilter mask the content of an input image according to the content of the input LabelMap . The masked pixel of the input image are set to the BackgroundValue. LabelMapMaskImageFilter can keep the input image for one label only, with Negated = false (the default) or it can mask the input image for a single label, when Negated equals true. In Both cases, the label is set with SetLabel() .

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


This implementation was taken from the Insight Journal paper: https://hdl.handle.net/1926/584 or http://www.insight-journal.org/browse/publication/176 

\see LabelMapToBinaryImageFilter , LabelMapToLabelImageFilter
\sa itk::simple::LabelMapMask for the procedural interface
\sa itk::LabelMapMaskImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT LabelMapMaskImageFilter : public ImageFilter {
    public:
      using Self = LabelMapMaskImageFilter;

      /** Destructor */
      virtual ~LabelMapMaskImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      LabelMapMaskImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = LabelPixelIDTypeList;

\

      /**
       * The label to mask or to not mask, depending on the value of the Negated ivar.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetLabel ( uint64_t Label ) { this->m_Label = Label; return *this; }

      /**
       * The label to mask or to not mask, depending on the value of the Negated ivar.
       */
      uint64_t GetLabel() const { return this->m_Label; }\

      /**
       * Set/Get the value used as "background" in the output image. Defaults to NumericTraits<PixelType>::ZeroValue() .
       */
      SITK_RETURN_SELF_TYPE_HEADER SetBackgroundValue ( double BackgroundValue ) { this->m_BackgroundValue = BackgroundValue; return *this; }

      /**
       * Set/Get the value used as "background" in the output image. Defaults to NumericTraits<PixelType>::ZeroValue() .
       */
      double GetBackgroundValue() const { return this->m_BackgroundValue; }\

      /**
       * Set/Get whether the Label should be masked or not.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetNegated ( bool Negated ) { this->m_Negated = Negated; return *this; }

      /** Set the value of Negated to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER NegatedOn() { return this->SetNegated(true); }
      SITK_RETURN_SELF_TYPE_HEADER NegatedOff() { return this->SetNegated(false); }

      /**
       * Set/Get whether the Label should be masked or not.
       */
      bool GetNegated() const { return this->m_Negated; }\

      /**
       * Set/Get whether the image size should be adjusted to the masked image or not.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetCrop ( bool Crop ) { this->m_Crop = Crop; return *this; }

      /** Set the value of Crop to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER CropOn() { return this->SetCrop(true); }
      SITK_RETURN_SELF_TYPE_HEADER CropOff() { return this->SetCrop(false); }

      /**
       * Set/Get whether the image size should be adjusted to the masked image or not.
       */
      bool GetCrop() const { return this->m_Crop; }\

      /**
       * Set/Get the boder added to the mask before the crop. The default is 0 on all the axes.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetCropBorder ( std::vector<unsigned int> CropBorder ) { this->m_CropBorder = std::move(CropBorder); return *this; }

      /** Set the values of the CropBorder vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetCropBorder( unsigned int value ) { this->m_CropBorder = std::vector<unsigned int>(3, value); return *this; }

      /**
       * Set/Get the boder added to the mask before the crop. The default is 0 on all the axes.
       */
      std::vector<unsigned int> GetCropBorder() const { return this->m_CropBorder; }

      /** Name of this class */
      std::string GetName() const { return std::string ("LabelMapMaskImageFilter"); }

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



      /*  */
      uint64_t  m_Label{1u};

      /*  */
      double  m_BackgroundValue{0};

      /*  */
      bool  m_Negated{false};

      /*  */
      bool  m_Crop{false};

      std::vector<unsigned int>  m_CropBorder{std::vector<unsigned int>(3, 0)};


    };

    /**\
     * \brief Mask and image with a LabelMap .
     *
     * This function directly calls the execute method of LabelMapMaskImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::LabelMapMaskImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image LabelMapMask ( const Image & labelMapImage, const Image & featureImage, uint64_t label = 1u, double backgroundValue = 0, bool negated = false, bool crop = false, std::vector<unsigned int> cropBorder = std::vector<unsigned int>(3, 0) );

     /** @} */
  }
}
#endif
