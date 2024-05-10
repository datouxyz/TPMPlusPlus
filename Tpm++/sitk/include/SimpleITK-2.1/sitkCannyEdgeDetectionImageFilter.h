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
#ifndef sitkCannyEdgeDetectionImageFilter_h
#define sitkCannyEdgeDetectionImageFilter_h

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

    /**\class CannyEdgeDetectionImageFilter
\brief This filter is an implementation of a Canny edge detector for scalar-valued images.

Based on John Canny's paper "A Computational Approach
 to Edge Detection"(IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. PAMI-8, No.6, November 1986), there are four major steps used in the edge-detection scheme: (1) Smooth the input image with Gaussian filter. (2) Calculate the second directional derivatives of the smoothed image. (3) Non-Maximum Suppression: the zero-crossings of 2nd derivative are found, and the sign of third derivative is used to find the correct extrema. (4) The hysteresis thresholding is applied to the gradient magnitude (multiplied with zero-crossings) of the smoothed image to find and link edges.

\par Inputs and Outputs
The input to this filter should be a scalar, real-valued Itk image of arbitrary dimension. The output should also be a scalar, real-value Itk image of the same dimensionality.


\par Parameters
There are four parameters for this filter that control the sub-filters used by the algorithm.


\par 
Variance and Maximum error are used in the Gaussian smoothing of the input image. See itkDiscreteGaussianImageFilter for information on these parameters.


\par 
Threshold is the lowest allowed value in the output image. Its data type is the same as the data type of the output image. Any values below the Threshold level will be replaced with the OutsideValue parameter value, whose default is zero.


TodoEdge-linking will be added when an itk connected component labeling algorithm is available.



\see DiscreteGaussianImageFilter 


\see ZeroCrossingImageFilter 


\see ThresholdImageFilter
\sa itk::simple::CannyEdgeDetection for the procedural interface
\sa itk::CannyEdgeDetectionImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT CannyEdgeDetectionImageFilter : public ImageFilter {
    public:
      using Self = CannyEdgeDetectionImageFilter;

      /** Destructor */
      virtual ~CannyEdgeDetectionImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      CannyEdgeDetectionImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = RealPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetLowerThreshold ( double LowerThreshold ) { this->m_LowerThreshold = LowerThreshold; return *this; }

      /**
       */
      double GetLowerThreshold() const { return this->m_LowerThreshold; }\

      /**
       * \brief Set the Threshold value for detected edges.
       * TODO: Document in the ITKv4 migration guide that the SetThreshold member function was removed from the CannyEdgeDetectionImageFilter , and that both UpperThreshold and LowerThreshold need to be set. To get the same results as with the SetThreshold method change "myfilter->SetThrehsold" to "myfilter->SetUpperThreshold", and add "myfilter->SetLowerThreshold(GetUpperThreshold()/2.0)"
       */
      SITK_RETURN_SELF_TYPE_HEADER SetUpperThreshold ( double UpperThreshold ) { this->m_UpperThreshold = UpperThreshold; return *this; }

      /**
       */
      double GetUpperThreshold() const { return this->m_UpperThreshold; }\

      /**
       * Set/Get the variance of the Gaussian smoothing filter.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetVariance ( std::vector<double> Variance ) { this->m_Variance = std::move(Variance); return *this; }

      /** Set the values of the Variance vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetVariance( double value ) { this->m_Variance = std::vector<double>(3, value); return *this; }

      /**
       * Set/Get the variance of the Gaussian smoothing filter.
       */
      std::vector<double> GetVariance() const { return this->m_Variance; }\

      /**
       * Set/Get the MaximumError parameter used by the Gaussian smoothing filter in this algorithm
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumError ( std::vector<double> MaximumError ) { this->m_MaximumError = std::move(MaximumError); return *this; }

      /** Set the values of the MaximumError vector all to value */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumError( double value ) { this->m_MaximumError = std::vector<double>(3, value); return *this; }

      /**
       * Set/Get the maximum error of the Gaussian smoothing kernel in each dimensional direction.
       */
      std::vector<double> GetMaximumError() const { return this->m_MaximumError; }

      /** Name of this class */
      std::string GetName() const { return std::string ("CannyEdgeDetectionImageFilter"); }

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


      double  m_LowerThreshold{0.0};

      double  m_UpperThreshold{0.0};

      /*  */
      std::vector<double>  m_Variance{std::vector<double>(3, 0.0)};

      /*  */
      std::vector<double>  m_MaximumError{std::vector<double>(3, 0.01)};


    };

    /**\
     * \brief This filter is an implementation of a Canny edge detector for scalar-valued images.
     *
     * This function directly calls the execute method of CannyEdgeDetectionImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::CannyEdgeDetectionImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image CannyEdgeDetection ( const Image& image1, double lowerThreshold = 0.0, double upperThreshold = 0.0, std::vector<double> variance = std::vector<double>(3, 0.0), std::vector<double> maximumError = std::vector<double>(3, 0.01) );

     /** @} */
  }
}
#endif