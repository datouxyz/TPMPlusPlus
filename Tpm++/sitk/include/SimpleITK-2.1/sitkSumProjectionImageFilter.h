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
#ifndef sitkSumProjectionImageFilter_h
#define sitkSumProjectionImageFilter_h

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

    /**\class SumProjectionImageFilter
\brief Sum projection.

This class was contributed to the Insight Journal by Gaetan Lehmann. The original paper can be found at https://hdl.handle.net/1926/164 

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


\see ProjectionImageFilter 


\see MedianProjectionImageFilter 


\see MeanProjectionImageFilter 


\see MeanProjectionImageFilter 


\see MaximumProjectionImageFilter 


\see MinimumProjectionImageFilter 


\see BinaryProjectionImageFilter 


\see StandardDeviationProjectionImageFilter
\sa itk::simple::SumProjection for the procedural interface
\sa itk::SumProjectionImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT SumProjectionImageFilter : public ImageFilter {
    public:
      using Self = SumProjectionImageFilter;

      /** Destructor */
      virtual ~SumProjectionImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      SumProjectionImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetProjectionDimension ( unsigned int ProjectionDimension ) { this->m_ProjectionDimension = ProjectionDimension; return *this; }

      /**
       */
      unsigned int GetProjectionDimension() const { return this->m_ProjectionDimension; }

      /** Name of this class */
      std::string GetName() const { return std::string ("SumProjectionImageFilter"); }

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


      unsigned int  m_ProjectionDimension{0u};


    };

    /**\
     * \brief Sum projection.
     *
     * This function directly calls the execute method of SumProjectionImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::SumProjectionImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image SumProjection ( const Image& image1, unsigned int projectionDimension = 0u );

     /** @} */
  }
}
#endif
