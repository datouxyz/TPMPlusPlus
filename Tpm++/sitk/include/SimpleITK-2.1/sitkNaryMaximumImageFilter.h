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
#ifndef sitkNaryMaximumImageFilter_h
#define sitkNaryMaximumImageFilter_h

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

   /**\class NaryMaximumImageFilter

\brief Computes the pixel-wise maximum of several images.

This class is templated over the types of the input images and the type of the output image. Numeric conversions (castings) are done by the C++ defaults.

The pixel type of the output images must have a valid definition of the operator<. This condition is required because internally this filter will perform an operation similar to:

\code
const OutputPixelType query_value = static_cast<OutputPixelType>(pixel_from_input_n);

if(current_maximum < query_value)

 {

 current_maximum = query_value;

 }

\endcode
 (where current_maximum is also of type OutputPixelType)

for each of the n input images.

For example, this filter could be used directly to find a "maximum projection" of a series of images, often used in preliminary analysis of time-series data.

\author Zachary Pincus


This filter was contributed by Zachary Pincus from the Department of Biochemistry and Program in Biomedical Informatics at Stanford University School of Medicine

\sa itk::simple::NaryMaximum for the procedural interface
   */
    class SITKBasicFilters_EXPORT NaryMaximumImageFilter
      : public ImageFilter
    {
    public:
      using Self = NaryMaximumImageFilter;

      /** Destructor */
      virtual ~NaryMaximumImageFilter();


      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      NaryMaximumImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = BasicPixelIDTypeList;




      /** Name of this class */
      std::string GetName() const { return std::string ("NaryMaximumImageFilter"); }

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



      bool m_InPlace{false};
    };


    /**
     * \brief Computes the pixel-wise maximum of several images.
     *
     * This function directly calls the execute method of NaryMaximumImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::NaryMaximumImageFilter for the object oriented interface
     * @{
     */
     SITKBasicFilters_EXPORT Image NaryMaximum ( const std::vector<Image> &images  );

     SITKBasicFilters_EXPORT Image NaryMaximum ( const Image& image1 );
     SITKBasicFilters_EXPORT Image NaryMaximum ( const Image& image1, const Image& image2 );
     SITKBasicFilters_EXPORT Image NaryMaximum ( const Image& image1, const Image& image2, const Image& image3 );
     SITKBasicFilters_EXPORT Image NaryMaximum ( const Image& image1, const Image& image2, const Image& image3, const Image& image4 );
     SITKBasicFilters_EXPORT Image NaryMaximum ( const Image& image1, const Image& image2, const Image& image3, const Image& image4, const Image& image5 );

     /** @{ */

}
}
#endif
