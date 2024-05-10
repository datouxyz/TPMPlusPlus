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
#ifndef sitkBoxSigmaImageFilter_h
#define sitkBoxSigmaImageFilter_h

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

    /**\class BoxSigmaImageFilter
\brief Implements a fast rectangular sigma filter using the accumulator approach.

This code was contributed in the Insight Journal paper: "Efficient implementation of kernel filtering" by Beare R., Lehmann G https://hdl.handle.net/1926/555 http://www.insight-journal.org/browse/publication/160 

\author Gaetan Lehmann
\sa itk::simple::BoxSigma for the procedural interface
\sa itk::BoxSigmaImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT BoxSigmaImageFilter : public ImageFilter {
    public:
      using Self = BoxSigmaImageFilter;

      /** Destructor */
      virtual ~BoxSigmaImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      BoxSigmaImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetRadius ( std::vector<unsigned int> Radius ) { this->m_Radius = std::move(Radius); return *this; }

      /** Set the values of the Radius vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetRadius( unsigned int value ) { this->m_Radius = std::vector<unsigned int>(3, value); return *this; }

      /**
       */
      std::vector<unsigned int> GetRadius() const { return this->m_Radius; }

      /** Name of this class */
      std::string GetName() const { return std::string ("BoxSigmaImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input image */

      Image Execute ( const Image& image1 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const Image& image1 );
      template <class TImageType> Image ExecuteInternal ( const Image& image1 );
      /** Dispatched methods which calls ExecuteInteral on each component */
      template <class TImageType> Image ExecuteInternalVectorImage ( const Image& image );

      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;
      friend struct detail::ExecuteInternalVectorImageAddressor<MemberFunctionType>;
      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      /*  */
      std::vector<unsigned int>  m_Radius{std::vector<unsigned int>(3, 1)};


    };

    /**\
     * \brief Implements a fast rectangular sigma filter using the accumulator approach.
     *
     * This function directly calls the execute method of BoxSigmaImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::BoxSigmaImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image BoxSigma ( const Image& image1, std::vector<unsigned int> radius = std::vector<unsigned int>(3, 1) );

     /** @} */
  }
}
#endif
