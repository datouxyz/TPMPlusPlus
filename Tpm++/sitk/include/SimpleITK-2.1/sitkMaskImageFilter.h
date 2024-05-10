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
#ifndef sitkMaskImageFilter_h
#define sitkMaskImageFilter_h

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

    /**\class MaskImageFilter
\brief Mask an image with a mask.

This class is templated over the types of the input image type, the mask image type and the type of the output image. Numeric conversions (castings) are done by the C++ defaults.

The pixel type of the input 2 image must have a valid definition of the operator != with zero. This condition is required because internally this filter will perform the operation

\code
if pixel_from_mask_image != masking_value

 pixel_output_image = pixel_input_image

else

 pixel_output_image = outside_value

\endcode


The pixel from the input 1 is cast to the pixel type of the output image.

Note that the input and the mask images must be of the same size.

\warning Any pixel value other than masking value (0 by default) will not be masked out.


\see MaskNegatedImageFilter
\sa itk::simple::Mask for the procedural interface
\sa itk::MaskImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT MaskImageFilter : public ImageFilter {
    public:
      using Self = MaskImageFilter;

      /** Destructor */
      virtual ~MaskImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      MaskImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = NonLabelPixelIDTypeList;

\

      /**
       * Method to explicitly set the outside value of the mask. Defaults to 0
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOutsideValue ( double OutsideValue ) { this->m_OutsideValue = OutsideValue; return *this; }

      /**
       */
      double GetOutsideValue() const { return this->m_OutsideValue; }\

      /**
       * Method to explicitly set the masking value of the mask. Defaults to 0
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaskingValue ( double MaskingValue ) { this->m_MaskingValue = MaskingValue; return *this; }

      /**
       * Method to get the masking value of the mask.
       */
      double GetMaskingValue() const { return this->m_MaskingValue; }

      /** Name of this class */
      std::string GetName() const { return std::string ("MaskImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */
#ifndef SWIG
      Image Execute ( Image && image, const Image & maskImage );
#endif
      Image Execute ( const Image & image, const Image & maskImage );


    private:
      /** Setup for member function dispatching */
      using MemberFunctionType = Image (Self::*)( const Image * image, const Image * maskImage );

      friend struct detail::DualExecuteInternalAddressor<MemberFunctionType>;
      template <class TImageType1, class TImageType2> Image DualExecuteInternal ( const Image * image, const Image * maskImage );


      std::unique_ptr<detail::DualMemberFunctionFactory<MemberFunctionType> > m_DualMemberFactory;



      double  m_OutsideValue{0};

      double  m_MaskingValue{0};


      bool m_InPlace{false};
    };

    /**\
     * \brief Mask an image with a mask.
     *
     * This function directly calls the execute method of MaskImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::MaskImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image Mask ( Image && image, const Image & maskImage, double outsideValue = 0, double maskingValue = 0 );
#endif
     SITKBasicFilters_EXPORT Image Mask ( const Image & image, const Image & maskImage, double outsideValue = 0, double maskingValue = 0 );

     /** @} */
  }
}
#endif