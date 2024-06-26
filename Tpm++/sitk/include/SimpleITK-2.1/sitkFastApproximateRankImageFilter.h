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
#ifndef sitkFastApproximateRankImageFilter_h
#define sitkFastApproximateRankImageFilter_h

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

    /**\class FastApproximateRankImageFilter
\brief A separable rank filter.

Medians aren't separable, but if you want a large robust smoother to be relatively quick then it is worthwhile pretending that they are.

This code was contributed in the Insight Journal paper: "Efficient implementation of kernel filtering" by Beare R., Lehmann G https://hdl.handle.net/1926/555 http://www.insight-journal.org/browse/publication/160 

\author Richard Beare
\sa itk::simple::FastApproximateRank for the procedural interface
\sa itk::FastApproximateRankImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT FastApproximateRankImageFilter : public ImageFilter {
    public:
      using Self = FastApproximateRankImageFilter;

      /** Destructor */
      virtual ~FastApproximateRankImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      FastApproximateRankImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetRank ( double Rank ) { this->m_Rank = Rank; return *this; }

      /**
       */
      double GetRank() const { return this->m_Rank; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetRadius ( std::vector<unsigned int> Radius ) { this->m_Radius = std::move(Radius); return *this; }

      /** Set the values of the Radius vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetRadius( unsigned int value ) { this->m_Radius = std::vector<unsigned int>(3, value); return *this; }

      /**
       */
      std::vector<unsigned int> GetRadius() const { return this->m_Radius; }

      /** Name of this class */
      std::string GetName() const { return std::string ("FastApproximateRankImageFilter"); }

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
      double  m_Rank{0.5};

      /*  */
      std::vector<unsigned int>  m_Radius{std::vector<unsigned int>(3, 1)};


    };

    /**\
     * \brief A separable rank filter.
     *
     * This function directly calls the execute method of FastApproximateRankImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::FastApproximateRankImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image FastApproximateRank ( const Image& image1, double rank = 0.5, std::vector<unsigned int> radius = std::vector<unsigned int>(3, 1) );

     /** @} */
  }
}
#endif
