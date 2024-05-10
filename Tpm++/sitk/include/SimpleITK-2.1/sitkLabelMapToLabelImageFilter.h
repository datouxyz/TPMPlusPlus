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
#ifndef sitkLabelMapToLabelImageFilter_h
#define sitkLabelMapToLabelImageFilter_h

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

    /**\class LabelMapToLabelImageFilter
\brief Converts a LabelMap to a labeled image.

LabelMapToBinaryImageFilter to a label image.

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


This implementation was taken from the Insight Journal paper: https://hdl.handle.net/1926/584 or http://www.insight-journal.org/browse/publication/176 

\see LabelMapToBinaryImageFilter , LabelMapMaskImageFilter
\sa itk::simple::LabelMapToLabel for the procedural interface
\sa itk::LabelMapToLabelImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT LabelMapToLabelImageFilter : public ImageFilter {
    public:
      using Self = LabelMapToLabelImageFilter;

      /** Destructor */
      virtual ~LabelMapToLabelImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      LabelMapToLabelImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = LabelPixelIDTypeList;


      /** Name of this class */
      std::string GetName() const { return std::string ("LabelMapToLabelImageFilter"); }

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



    };

    /**\
     * \brief Converts a LabelMap to a labeled image.
     *
     * This function directly calls the execute method of LabelMapToLabelImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::LabelMapToLabelImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image LabelMapToLabel ( const Image& image1 );

     /** @} */
  }
}
#endif
