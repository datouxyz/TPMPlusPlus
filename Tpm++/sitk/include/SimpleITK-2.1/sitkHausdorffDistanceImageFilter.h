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
#ifndef sitkHausdorffDistanceImageFilter_h
#define sitkHausdorffDistanceImageFilter_h

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

    /**\class HausdorffDistanceImageFilter
\brief Computes the Hausdorff distance between the set of non-zero pixels of two images.

HausdorffDistanceImageFilter computes the distance between the set non-zero pixels of two images using the following formula: \f[ H(A,B) = \max(h(A,B),h(B,A)) \f] where \f[ h(A,B) = \max_{a \in A} \min_{b \in B} \| a - b\| \f] is the directed Hausdorff distance and \f$A\f$ and \f$B\f$ are respectively the set of non-zero pixels in the first and second input images.

In particular, this filter uses the DirectedHausdorffImageFilter inside to compute the two directed distances and then select the largest of the two.

The Hausdorff distance measures the degree of mismatch between two sets and behaves like a metric over the set of all closed bounded sets - with properties of identity, symmetry and triangle inequality.

This filter requires the largest possible region of the first image and the same corresponding region in the second image. It behaves as filter with two inputs and one output. Thus it can be inserted in a pipeline with other filters. The filter passes the first input through unmodified.

This filter is templated over the two input image types. It assume both images have the same number of dimensions.

\see DirectedHausdorffDistanceImageFilter

\sa itk::HausdorffDistanceImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT HausdorffDistanceImageFilter : public ImageFilter {
    public:
      using Self = HausdorffDistanceImageFilter;

      /** Destructor */
      virtual ~HausdorffDistanceImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      HausdorffDistanceImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;

     /**
      * Return the computed Hausdorff distance.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetHausdorffDistance() const { return this->m_HausdorffDistance; };

     /**
      * Return the computed Hausdorff distance.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     double GetAverageHausdorffDistance() const { return this->m_AverageHausdorffDistance; };


      /** Name of this class */
      std::string GetName() const { return std::string ("HausdorffDistanceImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;


      /** Execute the filter on the input images */

      void Execute ( const Image& image1, const Image& image2 );

    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = void (Self::*)( const Image& image1, const Image& image2 );
      template <class TImageType> void ExecuteInternal ( const Image& image1, const Image& image2 );


      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      /* Some global documentation */
      double m_HausdorffDistance{0.0};
      /* Some global documentation */
      double m_AverageHausdorffDistance{0.0};


    };


  }
}
#endif
