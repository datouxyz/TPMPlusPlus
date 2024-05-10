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
#ifndef sitkScalarImageKmeansImageFilter_h
#define sitkScalarImageKmeansImageFilter_h

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

    /**\class ScalarImageKmeansImageFilter
\brief Classifies the intensity values of a scalar image using the K-Means algorithm.

Given an input image with scalar values, it uses the K-Means statistical classifier in order to define labels for every pixel in the image. The filter is templated over the type of the input image. The output image is predefined as having the same dimension of the input image and pixel type unsigned char, under the assumption that the classifier will generate less than 256 classes.

You may want to look also at the RelabelImageFilter that may be used as a postprocessing stage, in particular if you are interested in ordering the labels by their relative size in number of pixels.

\see Image 


\see ImageKmeansModelEstimator 


\see KdTreeBasedKmeansEstimator, WeightedCentroidKdTreeGenerator, KdTree 


\see RelabelImageFilter
\sa itk::simple::ScalarImageKmeans for the procedural interface
\sa itk::ScalarImageKmeansImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT ScalarImageKmeansImageFilter : public ImageFilter {
    public:
      using Self = ScalarImageKmeansImageFilter;

      /** Destructor */
      virtual ~ScalarImageKmeansImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      ScalarImageKmeansImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;
\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetClassWithInitialMean ( std::vector<double> ClassWithInitialMean ) { this->m_ClassWithInitialMean = ClassWithInitialMean; return *this; }

      /**
       */
      std::vector<double> GetClassWithInitialMean() const { return this->m_ClassWithInitialMean; }\

      /**
       * Set/Get the UseNonContiguousLabels flag. When this is set to false the labels are numbered contiguously, like in {0,1,3..N}. When the flag is set to true, the labels are selected in order to span the dynamic range of the output image. This last option is useful when the output image is intended only for display. The default value is false.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetUseNonContiguousLabels ( bool UseNonContiguousLabels ) { this->m_UseNonContiguousLabels = UseNonContiguousLabels; return *this; }

      /** Set the value of UseNonContiguousLabels to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER UseNonContiguousLabelsOn() { return this->SetUseNonContiguousLabels(true); }
      SITK_RETURN_SELF_TYPE_HEADER UseNonContiguousLabelsOff() { return this->SetUseNonContiguousLabels(false); }

      /**
       * Set/Get the UseNonContiguousLabels flag. When this is set to false the labels are numbered contiguously, like in {0,1,3..N}. When the flag is set to true, the labels are selected in order to span the dynamic range of the output image. This last option is useful when the output image is intended only for display. The default value is false.
       */
      bool GetUseNonContiguousLabels() const { return this->m_UseNonContiguousLabels; }
     /**
      * Return the array of Means found after the classification.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     std::vector<double> GetFinalMeans() const { return this->m_FinalMeans; };


      /** Name of this class */
      std::string GetName() const { return std::string ("ScalarImageKmeansImageFilter"); }

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


      std::vector<double>  m_ClassWithInitialMean{std::vector<double>()};

      /*  */
      bool  m_UseNonContiguousLabels{false};

      /* Docs */
      std::vector<double> m_FinalMeans{std::vector<double>()};


    };

    /**\
     * \brief Classifies the intensity values of a scalar image using the K-Means algorithm.
     *
     * This function directly calls the execute method of ScalarImageKmeansImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::ScalarImageKmeansImageFilter for the object oriented interface
     * @{
     */

     SITKBasicFilters_EXPORT Image ScalarImageKmeans ( const Image& image1, std::vector<double> classWithInitialMean = std::vector<double>(), bool useNonContiguousLabels = false );

     /** @} */
  }
}
#endif