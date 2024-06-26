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
#ifndef sitkRelabelComponentImageFilter_h
#define sitkRelabelComponentImageFilter_h

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

    /**\class RelabelComponentImageFilter
\brief Relabel the components in an image such that consecutive labels are used.

RelabelComponentImageFilter remaps the labels associated with the objects in an image (as from the output of ConnectedComponentImageFilter ) such that the label numbers are consecutive with no gaps between the label numbers used. By default, the relabeling will also sort the labels based on the size of the object: the largest object will have label #1, the second largest will have label #2, etc. If two labels have the same size their initial order is kept. The sorting by size can be disabled using SetSortByObjectSize.

Label #0 is assumed to be the background and is left unaltered by the relabeling.

RelabelComponentImageFilter is typically used on the output of the ConnectedComponentImageFilter for those applications that want to extract the largest object or the "k" largest objects. Any particular object can be extracted from the relabeled output using a BinaryThresholdImageFilter . A group of objects can be extracted from the relabled output using a ThresholdImageFilter .

Once all the objects are relabeled, the application can query the number of objects and the size of each object. Object sizes are returned in a vector. The size of the background is not calculated. So the size of object #1 is GetSizeOfObjectsInPixels() [0], the size of object #2 is GetSizeOfObjectsInPixels() [1], etc.

If user sets a minimum object size, all objects with fewer pixels than the minimum will be discarded, so that the number of objects reported will be only those remaining. The GetOriginalNumberOfObjects method can be called to find out how many objects were present before the small ones were discarded.

RelabelComponentImageFilter can be run as an "in place" filter, where it will overwrite its output. The default is run out of place (or generate a separate output). "In place" operation can be controlled via methods in the superclass, InPlaceImageFilter::InPlaceOn() and InPlaceImageFilter::InPlaceOff() .

\see ConnectedComponentImageFilter , BinaryThresholdImageFilter , ThresholdImageFilter
\sa itk::simple::RelabelComponent for the procedural interface
\sa itk::RelabelComponentImageFilter for the Doxygen on the original ITK class.
     */
    class SITKBasicFilters_EXPORT RelabelComponentImageFilter : public ImageFilter {
    public:
      using Self = RelabelComponentImageFilter;

      /** Destructor */
      virtual ~RelabelComponentImageFilter();

      /** Default Constructor that takes no arguments and initializes
       * default parameters */
      RelabelComponentImageFilter();

      /** Define the pixels types supported by this filter */
      using PixelIDTypeList = IntegerPixelIDTypeList;
\

      /**
       * Set the minimum size in pixels for an object. All objects smaller than this size will be discarded and will not appear in the output label map. NumberOfObjects will count only the objects whose pixel counts are greater than or equal to the minimum size. Call GetOriginalNumberOfObjects to find out how many objects were present in the original label map.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetMinimumObjectSize ( uint64_t MinimumObjectSize ) { this->m_MinimumObjectSize = MinimumObjectSize; return *this; }

      /**
       * Get the caller-defined minimum size of an object in pixels. If the caller has not set the minimum, 0 will be returned, which is to be interpreted as meaning that no minimum exists, and all objects in the original label map will be passed through to the output.
       */
      uint64_t GetMinimumObjectSize() const { return this->m_MinimumObjectSize; }\

      /**
       * Controls whether the object labels are sorted by size. If false, initial order of labels is kept.
       */
      SITK_RETURN_SELF_TYPE_HEADER SetSortByObjectSize ( bool SortByObjectSize ) { this->m_SortByObjectSize = SortByObjectSize; return *this; }

      /** Set the value of SortByObjectSize to true or false respectfully. */
      SITK_RETURN_SELF_TYPE_HEADER SortByObjectSizeOn() { return this->SetSortByObjectSize(true); }
      SITK_RETURN_SELF_TYPE_HEADER SortByObjectSizeOff() { return this->SetSortByObjectSize(false); }

      /**
       * Controls whether the object labels are sorted by size. If false, initial order of labels is kept.
       */
      bool GetSortByObjectSize() const { return this->m_SortByObjectSize; }
     /**
      * Get the number of objects in the image. This information is only valid after the filter has executed.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     uint32_t GetNumberOfObjects() const { return this->m_NumberOfObjects; };

     /**
      * Get the original number of objects in the image before small objects were discarded. This information is only valid after the filter has executed. If the caller has not specified a minimum object size, OriginalNumberOfObjects is the same as NumberOfObjects.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     uint32_t GetOriginalNumberOfObjects() const { return this->m_OriginalNumberOfObjects; };

     /**
      * Get the size of each object in physical space (in units of pixel size). This information is only valid after the filter has executed. Size of the background is not calculated. Size of object #1 is GetSizeOfObjectsInPhysicalUnits() [0]. Size of object #2 is GetSizeOfObjectsInPhysicalUnits() [1]. Etc.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     std::vector<float> GetSizeOfObjectsInPhysicalUnits() const { return this->m_SizeOfObjectsInPhysicalUnits; };

     /**
      * Get the size of each object in pixels. This information is only valid after the filter has executed. Size of the background is not calculated. Size of object #1 is GetSizeOfObjectsInPixels() [0]. Size of object #2 is GetSizeOfObjectsInPixels() [1]. Etc.
      *
      * This is a measurement. Its value is updated in the Execute
      * methods, so the value will only be valid after an execution.
      */
     std::vector<uint64_t> GetSizeOfObjectsInPixels() const { return this->m_SizeOfObjectsInPixels; };


      /** Name of this class */
      std::string GetName() const { return std::string ("RelabelComponentImageFilter"); }

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


      uint64_t  m_MinimumObjectSize{0u};

      bool  m_SortByObjectSize{true};

      /*  */
      uint32_t m_NumberOfObjects{0u};
      /*  */
      uint32_t m_OriginalNumberOfObjects{0u};
      /*  */
      std::vector<float> m_SizeOfObjectsInPhysicalUnits{std::vector<float>()};
      /*  */
      std::vector<uint64_t> m_SizeOfObjectsInPixels{std::vector<uint64_t>()};


      bool m_InPlace{false};
    };

    /**\
     * \brief Relabel the components in an image such that consecutive labels are used.
     *
     * This function directly calls the execute method of RelabelComponentImageFilter
     * in order to support a procedural API
     *
     * \sa itk::simple::RelabelComponentImageFilter for the object oriented interface
     * @{
     */
#ifndef SWIG
     SITKBasicFilters_EXPORT Image RelabelComponent ( Image&& image1, uint64_t minimumObjectSize = 0u, bool sortByObjectSize = true );
#endif
     SITKBasicFilters_EXPORT Image RelabelComponent ( const Image& image1, uint64_t minimumObjectSize = 0u, bool sortByObjectSize = true );

     /** @} */
  }
}
#endif
