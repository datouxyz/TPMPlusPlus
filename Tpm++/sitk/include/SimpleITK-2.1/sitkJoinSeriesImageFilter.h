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
#ifndef sitkJoinSeriesImageFilter_h
#define sitkJoinSeriesImageFilter_h

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

   /**\class JoinSeriesImageFilter

\brief Join N-D images into an (N+1)-D image.

This filter is templated over the input image type and the output image type. The pixel type of them must be the same and the input dimension must be less than the output dimension. When the input images are N-dimensional, they are joined in order and the size of the N+1'th dimension of the output is same as the number of the inputs. The spacing and the origin (where the first input is placed) for the N+1'th dimension is specified in this filter. The output image informations for the first N dimensions are taken from the first input. Note that all the inputs should have the same information.

\author Hideaki Hiraki


Contributed in the users list http://public.kitware.com/pipermail/insight-users/2004-February/006542.html

\sa itk::simple::JoinSeries for the procedural interface
   */
    class SITKBasicFilters_EXPORT JoinSeriesImageFilter
      : public ImageFilter
    {
    public:
      using Self = JoinSeriesImageFilter;

      /** Destructor */
      virtual ~JoinSeriesImageFilter();


      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      JoinSeriesImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = NonLabelPixelIDTypeList;


\

      /**
       * Set/Get origin of the new dimension
       */
      SITK_RETURN_SELF_TYPE_HEADER SetOrigin ( double Origin ) { this->m_Origin = Origin; return *this; }

      /**
       * Set/Get origin of the new dimension
       */
      double GetOrigin() const { return this->m_Origin; }\

      /**
       * Set/Get spacing of the new dimension
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSpacing ( double Spacing ) { this->m_Spacing = Spacing; return *this; }

      /**
       * Set/Get spacing of the new dimension
       */
      double GetSpacing() const { return this->m_Spacing; }

      /** Name of this class */
      std::string GetName() const { return std::string ("JoinSeriesImageFilter"); }

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


      double  m_Origin{0.0};

      double  m_Spacing{1.0};


    };


    /**
     * \brief Join N-D images into an (N+1)-D image.
     *
     * This function directly calls the execute method of JoinSeriesImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::JoinSeriesImageFilter for the object oriented interface
     * @{
     */
     SITKBasicFilters_EXPORT Image JoinSeries ( const std::vector<Image> &images , double origin = 0.0, double spacing = 1.0 );

     SITKBasicFilters_EXPORT Image JoinSeries ( const Image& image1, double origin = 0.0, double spacing = 1.0 );
     SITKBasicFilters_EXPORT Image JoinSeries ( const Image& image1, const Image& image2, double origin = 0.0, double spacing = 1.0 );
     SITKBasicFilters_EXPORT Image JoinSeries ( const Image& image1, const Image& image2, const Image& image3, double origin = 0.0, double spacing = 1.0 );
     SITKBasicFilters_EXPORT Image JoinSeries ( const Image& image1, const Image& image2, const Image& image3, const Image& image4, double origin = 0.0, double spacing = 1.0 );
     SITKBasicFilters_EXPORT Image JoinSeries ( const Image& image1, const Image& image2, const Image& image3, const Image& image4, const Image& image5, double origin = 0.0, double spacing = 1.0 );

     /** @{ */

}
}
#endif
