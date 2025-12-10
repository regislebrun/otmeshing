
#ifndef OTMESHING_PRIVATE_HXX
#define OTMESHING_PRIVATE_HXX

/* From http://gcc.gnu.org/wiki/Visibility */
/* Generic helper definitions for shared library support */
#if defined _WIN32 || defined __CYGWIN__
#define OTMESHING_HELPER_DLL_IMPORT __declspec(dllimport)
#define OTMESHING_HELPER_DLL_EXPORT __declspec(dllexport)
#define OTMESHING_HELPER_DLL_LOCAL
#else
#if __GNUC__ >= 4
#define OTMESHING_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
#define OTMESHING_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
#define OTMESHING_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define OTMESHING_HELPER_DLL_IMPORT
#define OTMESHING_HELPER_DLL_EXPORT
#define OTMESHING_HELPER_DLL_LOCAL
#endif
#endif

/* Now we use the generic helper definitions above to define OTMESHING_API and OTMESHING_LOCAL.
 * OTMESHING_API is used for the public API symbols. It either DLL imports or DLL exports (or does nothing for static build)
 * OTMESHING_LOCAL is used for non-api symbols. */

#ifndef OTMESHING_STATIC /* defined if OT is compiled as a DLL */
#ifdef OTMESHING_DLL_EXPORTS /* defined if we are building the OT DLL (instead of using it) */
#define OTMESHING_API OTMESHING_HELPER_DLL_EXPORT
#else
#define OTMESHING_API OTMESHING_HELPER_DLL_IMPORT
#endif /* OTMESHING_DLL_EXPORTS */
#define OTMESHING_LOCAL OTMESHING_HELPER_DLL_LOCAL
#else /* OTMESHING_STATIC is defined: this means OT is a static lib. */
#define OTMESHING_API
#define OTMESHING_LOCAL
#endif /* !OTMESHING_STATIC */


#endif // OTMESHING_PRIVATE_HXX

