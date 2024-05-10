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
#ifndef sitkVectorConfidenceConnectedImageFilter_h
#define sitkVectorConfidenceConnectedImageFilter_h

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

    /**\class VectorConfidenceConnectedImageFilter
\brief Segment pixels with similar statistics using connectivity.

This filter extracts a connected set of pixels whose pixel intensities are consistent with the pixel statistics of a seed point. The mean and variance across a neighborhood (8-connected, 26-connected, etc.) are calculated for a seed point. Then pixels connected to this seed point whose values are within the confidence interval for the seed point are grouped. The width of the confidence interval is controlled by the "Multiplier" variable (the confidence interval is the mean plus or minus the "Multiplier" times the standard deviation). If the intensity variations across a segment were gaussian, a "Multiplier" setting of 2.5 would define a confidence interval wide enough to capture 99% of samples in the segment.

After this initial segmentation is calculated, the mean and variance are re-calculated. All the pixels in the previous segmentation are used to calculate the mean the standard deviation (as opposed to using the pixels in the neighborhood of the seed point). The segmentation is then recalculted using these refined estimates for the mean and variance of the pixel values. This process is repeated for the specified number of iterations. Setting the "NumberOfIterations" to zero stops the algorithm after the initial segmentation from the seed point.

NOTE: the lower and upper threshold are restricted to lie within the valid numeric limits of the input data pixel type. Also, the limits may be adjusted to contain the seed point's intensity.
\sa itk::simple::VectorConfidenceConnected for the procedural interface
\sa itk::VectorConfidenceConnectedImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT VectorConfidenceConnectedImageFilter : public ImageFilter {
    public:
      using Self = VectorConfidenceConnectedImageFilter;

      /** Destructor */
      virtual ~VectorConfidenceConnectedImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      VectorConfidenceConnectedImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = VectorPixelIDTypeList;
\

      /**
       * \brief Set list of image indexes for seeds.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSeedList ( std::vector< std::vector<unsigned int> > SeedList ) { this->m_SeedList = std::move(SeedList); return *this; }

      /**
       * \brief Get list of seeds
       */
      std::vector< std::vector< unsigned int > > GetSeedList() const { return this->m_SeedList; }
      /** \brief Add SeedList point */
      SITK_RETURN_SELF_TYPE_HEADER AddSeed( std::vector< unsigned int > point ) { this->m_SeedList.push_back(std::move(point)); return *this;}
      /** \brief Remove all SeedList points */
      SITK_RETURN_SELF_TYPE_HEADER ClearSeeds( ) { this->m_SeedList.clear(); return *this;}
\

      /**
       * Set/Get the number of iterations
       */
      SITK_RETURN_SELF_TYPE_HEADER SetNumberOfIterations ( unsigned int NumberOfIterations ) { this->m_NumberOfIterations = NumberOfIterations; return *this; }

      /**
       * Set/Get the number of iterations
       */
      unsigned int GetNumberOfIterations() const { return this->m_NumberOfIterations; }\

      /**
       * Set/Get the multiplier to define the confidence interval. Multiplier can be anything greater than zero. A typical value is 2.5
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMultiplier ( double Multiplier ) { this->m_Multiplier = Multiplier; return *this; }

      /**
       * Set/Get the multiplier to define the confidence interval. Multiplier can be anything greater than zero. A typical value is 2.5
       */
      double GetMultiplier() const { return this->m_Multiplier; }\

      /**
       * Get/Set the radius of the neighborhood over which the statistics are evaluated
       */
      SITK_RETURN_SELF_TYPE_HEADER SetInitialNeighborhoodRadius ( unsigned int InitialNeighborhoodRadius ) { this->m_InitialNeighborhoodRadius = InitialNeighborhoodRadius; return *this; }

      /**
       * Get/Set the radius of the neighborhood over which the statistics are evaluated
       */
      unsigned int GetInitialNeighborhoodRadius() const { return this->m_InitialNeighborhoodRadius; }\

      /**
       * Set/Get value to replace thresholded pixels
       */
      SITK_RETURN_SELF_TYPE_HEADER SetReplaceValue ( uint8_t ReplaceValue ) { this->m_ReplaceValue = ReplaceValue; return *this; }

      /**
       * Set/Get value to replace thresholded pixels
       */
      uint8_t GetReplaceValue() const { return this->m_ReplaceValue; }
     /**
      * Get the Mean Vector computed during the segmentation
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     std::vector<double> GetMean() const { return this->m_Mean; };

     /**
      * Get the Covariance matrix computed during the segmentation
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     std::vector<double> GetCovariance() const { return this->m_Covariance; };


      /** Name of this class */
      std::string GetName() const { return std::string ("VectorConfidenceConnectedImageFilter"); }

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


      std::vector< std::vector<unsigned int> >  m_SeedList{std::vector< std::vector<unsigned int > >()};

      unsigned int  m_NumberOfIterations{4u};

      double  m_Multiplier{4.5};

      unsigned int  m_InitialNeighborhoodRadius{1u};

      uint8_t  m_ReplaceValue{1u};

      /* Some global documentation */
      std::vector<double> m_Mean{std::vector<double>()};
      /* Some global documentation */
      std::vector<double> m_Covariance{std::vector<double>()};


    };

    /**\
     * \brief Segment pixels with similar statistics using connectivity.
     *
     * This function directly calls the execute method of VectorConfidenceConnectedImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::VectorConfidenceConnectedImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image VectorConfidenceConnected ( const Image& image1, std::vector< std::vector<unsigned int> > seedList = std::vector< std::vector<unsigned int > >(), unsigned int numberOfIterations = 4u, double multiplier = 4.5, unsigned int initialNeighborhoodRadius = 1u, uint8_t replaceValue = 1u );

     /** @} */
  }
}
#endif