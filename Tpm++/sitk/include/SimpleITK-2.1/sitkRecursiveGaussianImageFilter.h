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
#ifndef sitkRecursiveGaussianImageFilter_h
#define sitkRecursiveGaussianImageFilter_h

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

    /**\class RecursiveGaussianImageFilter
\brief Base class for computing IIR convolution with an approximation of a Gaussian kernel.

\f[ \frac{ 1 }{ \sigma \sqrt{ 2 \pi } } \exp{ \left( - \frac{x^2}{ 2 \sigma^2 } \right) } \f] 

RecursiveGaussianImageFilter is the base class for recursive filters that approximate convolution with the Gaussian kernel. This class implements the recursive filtering method proposed by R.Deriche in IEEE-PAMI Vol.12, No.1, January 1990, pp 78-87, "Fast Algorithms for Low-Level Vision"

Details of the implementation are described in the technical report: R. Deriche, "Recursively Implementing The Gaussian and Its Derivatives", INRIA, 1993, ftp://ftp.inria.fr/INRIA/tech-reports/RR/RR-1893.ps.gz 

Further improvements of the algorithm are described in: G. Farneback & C.-F. Westin, "On Implementation of Recursive Gaussian
 Filters", so far unpublished.

As compared to itk::DiscreteGaussianImageFilter , this filter tends to be faster for large kernels, and it can take the derivative of the blurred image in one step. Also, note that we have itk::RecursiveGaussianImageFilter::SetSigma() , but itk::DiscreteGaussianImageFilter::SetVariance() .

\see DiscreteGaussianImageFilter
\sa itk::simple::RecursiveGaussian for the procedural interface
\sa itk::RecursiveGaussianImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT RecursiveGaussianImageFilter : public ImageFilter {
    public:
      using Self = RecursiveGaussianImageFilter;

      /** Destructor */
      virtual ~RecursiveGaussianImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      RecursiveGaussianImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = typelist::Append<BasicPixelIDTypeList, VectorPixelIDTypeList>::Type;
\

      /**
       * Set/Get the Sigma, measured in world coordinates, of the Gaussian kernel. The default is 1.0. An exception will be generated if the Sigma value is less than or equal to zero.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSigma ( double Sigma ) { this->m_Sigma = Sigma; return *this; }

      /**
       * Set/Get the Sigma, measured in world coordinates, of the Gaussian kernel. The default is 1.0. An exception will be generated if the Sigma value is less than or equal to zero.
       */
      double GetSigma() const { return this->m_Sigma; }\

      /**
       * Set/Get the flag for normalizing the gaussian over scale-space.

This flag enables the analysis of the differential shape of features independent of their size ( both pixels and physical size ). Following the notation of Tony Lindeberg:

Let \f[ L(x; t) = g(x; t) \ast f(x) \f] be the scale-space representation of image \f[ f(x) \f] where \f[ g(x; t) = \frac{1}{ \sqrt{ 2 \pi t} } \exp{ \left( -\frac{x^2}{ 2 t } \right) } \f] is the Gaussian function and \f[\ast\f] denotes convolution. This is a change from above with \f[ t = \sigma^2 \f] .

Then the normalized derivative operator for normalized coordinates across scale is:

 \f[ \partial_\xi = \sqrt{t} \partial_x \f] 

The resulting scaling factor is \f[ \sigma^N \f] where N is the order of the derivative.

When this flag is ON the filter will be normalized in such a way that the values of derivatives are not biased by the size of the object. That is to say the maximum value a feature reaches across scale is independent of the scale of the object.

For analyzing an image across scale-space you want to enable this flag. It is disabled by default.

\note Not all scale space axioms are satisfied by this filter, some are only approximated. Particularly, at fine scales ( say less than 1 pixel ) other methods such as a discrete Gaussian kernel should be considered.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetNormalizeAcrossScale ( bool NormalizeAcrossScale ) { this->m_NormalizeAcrossScale = NormalizeAcrossScale; return *this; }

      /** Set the value of NormalizeAcrossScale to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER NormalizeAcrossScaleOn() { return this->SetNormalizeAcrossScale(true); }
      SITK_RETURN_SELF_TYPE_HEADER NormalizeAcrossScaleOff() { return this->SetNormalizeAcrossScale(false); }

      /**
       */
      bool GetNormalizeAcrossScale() const { return this->m_NormalizeAcrossScale; }

      typedef enum {ZeroOrder,FirstOrder,SecondOrder} OrderType;\

      /**
       * Set/Get the Order of the Gaussian to convolve with. 

\li ZeroOrder is equivalent to convolving with a Gaussian. This is the default. 

\li FirstOrder is equivalent to convolving with the first derivative of a Gaussian. 

\li SecondOrder is equivalent to convolving with the second derivative of a Gaussian.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOrder ( OrderType Order ) { this->m_Order = Order; return *this; }

      /**
       * Set/Get the Order of the Gaussian to convolve with. 

\li ZeroOrder is equivalent to convolving with a Gaussian. This is the default. 

\li FirstOrder is equivalent to convolving with the first derivative of a Gaussian. 

\li SecondOrder is equivalent to convolving with the second derivative of a Gaussian.
       */
      OrderType GetOrder() const { return this->m_Order; }\

      /**
       */
      SITK_RETURN_SELF_TYPE_HEADER SetDirection ( unsigned int Direction ) { this->m_Direction = Direction; return *this; }

      /**
       */
      unsigned int GetDirection() const { return this->m_Direction; }

      /** Name of this class */
      std::string GetName() const { return std::string ("RecursiveGaussianImageFilter"); }

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


      /*  */
      double  m_Sigma{1.0};

      /*  */
      bool  m_NormalizeAcrossScale{false};

      /*  */
      OrderType  m_Order{itk::simple::RecursiveGaussianImageFilter::ZeroOrder};

      /*  */
      unsigned int  m_Direction{0u};


      bool m_InPlace{false};
    };

    /**\
     * \brief Base class for computing IIR convolution with an approximation of a Gaussian kernel.
     *
     * This function directly calls the execute method of RecursiveGaussianImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::RecursiveGaussianImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image RecursiveGaussian ( Image&& image1, double sigma = 1.0, bool normalizeAcrossScale = false, RecursiveGaussianImageFilter::OrderType order = itk::simple::RecursiveGaussianImageFilter::ZeroOrder, unsigned int direction = 0u );
#endif
     SITKBasicFilters_EXPORT Image RecursiveGaussian ( const Image& image1, double sigma = 1.0, bool normalizeAcrossScale = false, RecursiveGaussianImageFilter::OrderType order = itk::simple::RecursiveGaussianImageFilter::ZeroOrder, unsigned int direction = 0u );

     /** @} */
  }
}
#endif
