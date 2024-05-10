
#ifndef ITKIOTransformMINC_EXPORT_H
#define ITKIOTransformMINC_EXPORT_H

#ifdef ITK_STATIC
#  define ITKIOTransformMINC_EXPORT
#  define ITKIOTransformMINC_HIDDEN
#else
#  ifndef ITKIOTransformMINC_EXPORT
#    ifdef ITKIOTransformMINC_EXPORTS
        /* We are building this library */
#      define ITKIOTransformMINC_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define ITKIOTransformMINC_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef ITKIOTransformMINC_HIDDEN
#    define ITKIOTransformMINC_HIDDEN 
#  endif
#endif

#ifndef ITKIOTRANSFORMMINC_DEPRECATED
#  define ITKIOTRANSFORMMINC_DEPRECATED __declspec(deprecated)
#endif

#ifndef ITKIOTRANSFORMMINC_DEPRECATED_EXPORT
#  define ITKIOTRANSFORMMINC_DEPRECATED_EXPORT ITKIOTransformMINC_EXPORT ITKIOTRANSFORMMINC_DEPRECATED
#endif

#ifndef ITKIOTRANSFORMMINC_DEPRECATED_NO_EXPORT
#  define ITKIOTRANSFORMMINC_DEPRECATED_NO_EXPORT ITKIOTransformMINC_HIDDEN ITKIOTRANSFORMMINC_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef ITKIOTRANSFORMMINC_NO_DEPRECATED
#    define ITKIOTRANSFORMMINC_NO_DEPRECATED
#  endif
#endif

#endif /* ITKIOTransformMINC_EXPORT_H */
