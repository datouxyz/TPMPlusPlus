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
#ifndef sitkMergeLabelMapFilter_h
#define sitkMergeLabelMapFilter_h

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

   /**\class MergeLabelMapFilter

\brief Merges several Label Maps.

This filter takes one or more input Label Map and merges them.

SetMethod() can be used to change how the filter manage the labels from the different label maps. KEEP (0): MergeLabelMapFilter do its best to keep the label unchanged, but if a label is already used in a previous label map, a new label is assigned. AGGREGATE (1): If the same label is found several times in the label maps, the label objects with the same label are merged. PACK (2): MergeLabelMapFilter relabel all the label objects by order of processing. No conflict can occur. STRICT (3): MergeLabelMapFilter keeps the labels unchanged and raises an exception if the same label is found in several images.

This implementation was taken from the Insight Journal paper: https://hdl.handle.net/1926/584 or http://www.insight-journal.org/browse/publication/176 

\author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.


\see ShapeLabelObject , RelabelComponentImageFilter

\sa itk::simple::MergeLabelMapFilter for the procedural interface
   */
    class SITKBasicFilters_EXPORT MergeLabelMapFilter
      : public ImageFilter
    {
    public:
      using Self = MergeLabelMapFilter;

      /** Destructor */
      virtual ~MergeLabelMapFilter();


      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      MergeLabelMapFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = LabelPixelIDTypeList;




      typedef enum {Keep,Aggregate,Pack,Strict} MethodType;\

      /**
       * Set/Get the method used to merge the label maps
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMethod ( MethodType Method ) { this->m_Method = Method; return *this; }

      /**
       * Set/Get the method used to merge the label maps
       */
      MethodType GetMethod() const { return this->m_Method; }

      /** Name of this class */
      std::string GetName() const { return std::string ("MergeLabelMapFilter"); }

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


      MethodType  m_Method{itk::simple::MergeLabelMapFilter::Keep};


      bool m_InPlace{false};
    };


    /**
     * \brief Merges several Label Maps.
     *
     * This function directly calls the execute method of MergeLabelMapFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::MergeLabelMapFilter for the object oriented interface
     * @{
     */
     SITKBasicFilters_EXPORT Image MergeLabelMap ( const std::vector<Image> &images , MergeLabelMapFilter::MethodType method = itk::simple::MergeLabelMapFilter::Keep );

     SITKBasicFilters_EXPORT Image MergeLabelMap ( const Image& image1, MergeLabelMapFilter::MethodType method = itk::simple::MergeLabelMapFilter::Keep );
     SITKBasicFilters_EXPORT Image MergeLabelMap ( const Image& image1, const Image& image2, MergeLabelMapFilter::MethodType method = itk::simple::MergeLabelMapFilter::Keep );
     SITKBasicFilters_EXPORT Image MergeLabelMap ( const Image& image1, const Image& image2, const Image& image3, MergeLabelMapFilter::MethodType method = itk::simple::MergeLabelMapFilter::Keep );
     SITKBasicFilters_EXPORT Image MergeLabelMap ( const Image& image1, const Image& image2, const Image& image3, const Image& image4, MergeLabelMapFilter::MethodType method = itk::simple::MergeLabelMapFilter::Keep );
     SITKBasicFilters_EXPORT Image MergeLabelMap ( const Image& image1, const Image& image2, const Image& image3, const Image& image4, const Image& image5, MergeLabelMapFilter::MethodType method = itk::simple::MergeLabelMapFilter::Keep );

     /** @{ */

}
}
#endif
