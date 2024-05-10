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
#ifndef sitkChangeLabelImageFilter_h
#define sitkChangeLabelImageFilter_h

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

    /**\class ChangeLabelImageFilter
\brief Change Sets of Labels.

This filter produces an output image whose pixels are either copied from the input if they are not being changed or are rewritten based on the change parameters

This filter is templated over the input image type and the output image type.

The filter expect both images to have the same number of dimensions.

\author Tim Kelliher. GE Research, Niskayuna, NY. 


\note This work was supported by a grant from DARPA, executed by the U.S. Army Medical Research and Materiel Command/TATRC Assistance Agreement, Contract::W81XWH-05-2-0059.
\sa itk::simple::ChangeLabel for the procedural interface
\sa itk::ChangeLabelImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT ChangeLabelImageFilter : public ImageFilter {
    public:
      using Self = ChangeLabelImageFilter;

      /** Destructor */
      virtual ~ChangeLabelImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      ChangeLabelImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;
\

      /**
       * Set the entire change map
       */
      SITK_RETURN_SELF_TYPE_HEADER SetChangeMap ( std::map<double,double> ChangeMap ) { this->m_ChangeMap = ChangeMap; return *this; }

      /**
       */
      std::map<double,double> GetChangeMap() const { return this->m_ChangeMap; }

      /** Name of this class */
      std::string GetName() const { return std::string ("ChangeLabelImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */
#ifndef SWIG
      Image Execute ( Image&& image1 );
#endif
      Image Execute ( const Image& image1 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      std::map<double,double>  m_ChangeMap{std::map<double,double>()};


      bool m_InPlace{false};
    };

    /**\
     * \brief Change Sets of Labels.
     *
     * This function directly calls the execute method of ChangeLabelImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::ChangeLabelImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image ChangeLabel ( Image&& image1, std::map<double,double> changeMap = std::map<double,double>() );
#endif
     SITKBasicFilters_EXPORT Image ChangeLabel ( const Image& image1, std::map<double,double> changeMap = std::map<double,double>() );

     /** @} */
  }
}
#endif