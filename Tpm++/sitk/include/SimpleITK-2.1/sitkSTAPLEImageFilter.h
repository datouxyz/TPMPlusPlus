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
#ifndef sitkSTAPLEImageFilter_h
#define sitkSTAPLEImageFilter_h

/*
 * WARNING: DO NOT EDIT THIS FILE!
 * THIS FILE IS AUTOMATICALLY GENERATED BY THE SIMPLEITK BUILD PROCESS.
 * Please look at sitkMultiInputImageFilterTemplate.h.in to make changes.
 */

#include <memory>

#include "sitkBasicFilters.h"
#include "sitkImageFilter.h"

namespace itk {
  namespace simple {

   /**\class STAPLEImageFilter

\brief The STAPLE filter implements the Simultaneous Truth and Performance Level Estimation algorithm for generating ground truth volumes from a set of binary expert segmentations.

The STAPLE algorithm treats segmentation as a pixelwise classification, which leads to an averaging scheme that accounts for systematic biases in the behavior of experts in order to generate a fuzzy ground truth volume and simultaneous accuracy assessment of each expert. The ground truth volumes produced by this filter are floating point volumes of values between zero and one that indicate probability of each pixel being in the object targeted by the segmentation.

The STAPLE algorithm is described in

S. Warfield, K. Zou, W. Wells, "Validation of image segmentation and expert
 quality with an expectation-maximization algorithm" in MICCAI 2002: Fifth International Conference on Medical Image Computing and Computer-Assisted Intervention, Springer-Verlag, Heidelberg, Germany, 2002, pp. 298-306

\par INPUTS
Input volumes to the STAPLE filter must be binary segmentations of an image, that is, there must be a single foreground value that represents positively classified pixels (pixels that are considered to belong inside the segmentation). Any number of background pixel values may be present in the input images. You can, for example, input volumes with many different labels as long as the structure you are interested in creating ground truth for is consistently labeled among all input volumes. Pixel type of the input volumes does not matter. Specify the label value for positively classified pixels using SetForegroundValue. All other labels will be considered to be negatively classified pixels (background).


Input volumes must all contain the same size RequestedRegions.

\par OUTPUTS
The STAPLE filter produces a single output volume with a range of floating point values from zero to one. IT IS VERY IMPORTANT TO INSTANTIATE THIS FILTER WITH A FLOATING POINT OUTPUT TYPE (floats or doubles). You may threshold the output above some probability threshold if you wish to produce a binary ground truth.


\par PARAMETERS
The STAPLE algorithm requires a number of inputs. You may specify any number of input volumes using the SetInput(i, p_i) method, where i ranges from zero to N-1, N is the total number of input segmentations, and p_i is the SmartPointer to the i-th segmentation.


The SetConfidenceWeight parameter is a modifier for the prior probability that any pixel would be classified as inside the target object. This implementation of the STAPLE algorithm automatically calculates prior positive classification probability as the average fraction of the image volume filled by the target object in each input segmentation. The ConfidenceWeight parameter allows for scaling the of this default prior probability: if g_t is the prior probability that a pixel would be classified inside the target object, then g_t is set to g_t * ConfidenceWeight before iterating on the solution. In general ConfidenceWeight should be left to the default of 1.0.

You must provide a foreground value using SetForegroundValue that the STAPLE algorithm will use to identify positively classified pixels in the the input images. All other values in the image will be treated as background values. For example, if your input segmentations consist of 1's everywhere inside the segmented region, then use SetForegroundValue(1).

The STAPLE algorithm is an iterative E-M algorithm and will converge on a solution after some number of iterations that cannot be known a priori. After updating the filter, the total elapsed iterations taken to converge on the solution can be queried through GetElapsedIterations() . You may also specify a MaximumNumberOfIterations, after which the algorithm will stop iterating regardless of whether or not it has converged. This implementation of the STAPLE algorithm will find the solution to within seven digits of precision unless it is stopped early.

Once updated, the Sensitivity (true positive fraction, q) and Specificity (true negative fraction, q) for each expert input volume can be queried using GetSensitivity(i) and GetSpecificity(i), where i is the i-th input volume.

\par REQUIRED PARAMETERS
The only required parameters for this filter are the ForegroundValue and the input volumes. All other parameters may be safely left to their default values. Please see the paper cited above for more information on the STAPLE algorithm and its parameters. A proper understanding of the algorithm is important for interpreting the results that it produces.


\par EVENTS
This filter invokes IterationEvent() at each iteration of the E-M algorithm. Setting the AbortGenerateData() flag will cause the algorithm to halt after the current iteration and produce results just as if it had converged. The algorithm makes no attempt to report its progress since the number of iterations needed cannot be known in advance.

\sa itk::simple::STAPLE for the procedural interface
   */
    class SITKBasicFilters_EXPORT STAPLEImageFilter
      : public ImageFilter
    {
    public:
      using Self = STAPLEImageFilter;

      /** Destructor */
      virtual ~STAPLEImageFilter();


      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      STAPLEImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;


\

      /**
       * Scales the estimated prior probability that a pixel will be inside the targeted object of segmentation. The default prior probability g_t is calculated automatically as the average fraction of positively classified pixels to the total size of the volume (across all input volumes). ConfidenceWeight will scale this default value as g_t = g_t * ConfidenceWeight. In general, ConfidenceWeight should be left to the default of 1.0.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetConfidenceWeight ( double ConfidenceWeight ) { this->m_ConfidenceWeight = ConfidenceWeight; return *this; }

      /**
       * Scales the estimated prior probability that a pixel will be inside the targeted object of segmentation. The default prior probability g_t is calculated automatically as the average fraction of positively classified pixels to the total size of the volume (across all input volumes). ConfidenceWeight will scale this default value as g_t = g_t * ConfidenceWeight. In general, ConfidenceWeight should be left to the default of 1.0.
       */
      double GetConfidenceWeight() const { return this->m_ConfidenceWeight; }\

      /**
       * Set get the binary ON value of the input image.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetForegroundValue ( double ForegroundValue ) { this->m_ForegroundValue = ForegroundValue; return *this; }

      /**
       * Set get the binary ON value of the input image.
       */
      double GetForegroundValue() const { return this->m_ForegroundValue; }\

      /**
       * Set/Get the maximum number of iterations after which the STAPLE algorithm will be considered to have converged. In general this SHOULD NOT be set and the algorithm should be allowed to converge on its own.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMaximumIterations ( unsigned int MaximumIterations ) { this->m_MaximumIterations = MaximumIterations; return *this; }

      /**
       * Set/Get the maximum number of iterations after which the STAPLE algorithm will be considered to have converged. In general this SHOULD NOT be set and the algorithm should be allowed to converge on its own.
       */
      unsigned int GetMaximumIterations() const { return this->m_MaximumIterations; }
     /**
      * Get the number of elapsed iterations of the iterative E-M algorithm.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     uint32_t GetElapsedIterations() const { return this->m_ElapsedIterations; };

     /**
      * After the filter is updated, this method returns a std::vector<double> of all Sensitivity (true positive fraction, p) values for the expert input volumes.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     std::vector<double> GetSensitivity() const { return this->m_Sensitivity; };

     /** \brief After the filter is updated, this method returns the Specificity (true negative fraction, q) value for the i-th expert input volume.
      *
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     std::vector<double> GetSpecificity () const { return this->m_Specificity ; };


      /** Name of this class */
      std::string GetName() const { return std::string ("STAPLEImageFilter"); }

      /** Print ourselves out */
      std::string ToString() const;

      /** Execute the filter on the input images */
      Image Execute ( const std::vector<Image> &images);
      Image Execute ( const Image& image1 );
      Image Execute ( const Image& image1, const Image& image2 );
      Image Execute ( const Image& image1, const Image& image2, const Image& image3 );
      Image Execute ( const Image& image1, const Image& image2, const Image& image3, const Image& image4 );
      Image Execute ( const Image& image1, const Image& image2, const Image& image3, const Image& image4, const Image& image5 );




    private:

      /** Setup for member function dispatching */

      using MemberFunctionType = Image (Self::*)( const std::vector<Image> & );
      template <class TImageType> Image ExecuteInternal ( const std::vector<Image> &images );



      friend struct detail::MemberFunctionAddressor<MemberFunctionType>;

      std::unique_ptr<detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


      double  m_ConfidenceWeight{1.0};

      double  m_ForegroundValue{1.0};

      unsigned int  m_MaximumIterations{std::numeric_limits<unsigned int>::max()};


      uint32_t m_ElapsedIterations{0};

      std::vector<double> m_Sensitivity{std::vector<double>()};

      std::vector<double> m_Specificity {std::vector<double>()};


    };


    /**
     * \brief The STAPLE filter implements the Simultaneous Truth and Performance Level Estimation algorithm for generating ground truth volumes from a set of binary expert segmentations.
     *
     * This function directly calls the execute method of STAPLEImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::STAPLEImageFilter for the object oriented interface
     * @{
     */
     SITKBasicFilters_EXPORT Image STAPLE ( const std::vector<Image> &images , double confidenceWeight = 1.0, double foregroundValue = 1.0, unsigned int maximumIterations = std::numeric_limits<unsigned int>::max() );

     SITKBasicFilters_EXPORT Image STAPLE ( const Image& image1, double confidenceWeight = 1.0, double foregroundValue = 1.0, unsigned int maximumIterations = std::numeric_limits<unsigned int>::max() );
     SITKBasicFilters_EXPORT Image STAPLE ( const Image& image1, const Image& image2, double confidenceWeight = 1.0, double foregroundValue = 1.0, unsigned int maximumIterations = std::numeric_limits<unsigned int>::max() );
     SITKBasicFilters_EXPORT Image STAPLE ( const Image& image1, const Image& image2, const Image& image3, double confidenceWeight = 1.0, double foregroundValue = 1.0, unsigned int maximumIterations = std::numeric_limits<unsigned int>::max() );
     SITKBasicFilters_EXPORT Image STAPLE ( const Image& image1, const Image& image2, const Image& image3, const Image& image4, double confidenceWeight = 1.0, double foregroundValue = 1.0, unsigned int maximumIterations = std::numeric_limits<unsigned int>::max() );
     SITKBasicFilters_EXPORT Image STAPLE ( const Image& image1, const Image& image2, const Image& image3, const Image& image4, const Image& image5, double confidenceWeight = 1.0, double foregroundValue = 1.0, unsigned int maximumIterations = std::numeric_limits<unsigned int>::max() );

     /** @{ */

}
}
#endif
