
#ifndef SimpleITKFilters_EXPORT_H
#define SimpleITKFilters_EXPORT_H

#ifdef ITK_STATIC
#  define SimpleITKFilters_EXPORT
#  define SimpleITKFilters_HIDDEN
#else
#  ifndef SimpleITKFilters_EXPORT
#    ifdef SimpleITKFilters_EXPORTS
        /* We are building this library */
#      define SimpleITKFilters_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define SimpleITKFilters_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef SimpleITKFilters_HIDDEN
#    define SimpleITKFilters_HIDDEN 
#  endif
#endif

#ifndef SIMPLEITKFILTERS_DEPRECATED
#  define SIMPLEITKFILTERS_DEPRECATED __declspec(deprecated)
#endif

#ifndef SIMPLEITKFILTERS_DEPRECATED_EXPORT
#  define SIMPLEITKFILTERS_DEPRECATED_EXPORT SimpleITKFilters_EXPORT SIMPLEITKFILTERS_DEPRECATED
#endif

#ifndef SIMPLEITKFILTERS_DEPRECATED_NO_EXPORT
#  define SIMPLEITKFILTERS_DEPRECATED_NO_EXPORT SimpleITKFilters_HIDDEN SIMPLEITKFILTERS_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef SIMPLEITKFILTERS_NO_DEPRECATED
#    define SIMPLEITKFILTERS_NO_DEPRECATED
#  endif
#endif

#endif /* SimpleITKFilters_EXPORT_H */
