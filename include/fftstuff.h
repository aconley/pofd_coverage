#ifndef __fftstuff__
#define __fftstuff__

#include<fftw3.h>

// FFT precision type
#if USEFFT32BIT

#define FFTFLOAT float
#define FITSIOIMG FLOAT_IMG
#define TFITSFLOAT TFLOAT
#define FFTWMALLOC fftwf_malloc
#define FFTWFREE fftwf_free
#define FFTWCOMPLEX fftwf_complex
#define FFTWPLAN fftwf_plan
#define FFTWDFTR2C1D fftwf_plan_dft_r2c_1d
#define FFTWDFTC2R1D fftwf_plan_dft_c2r_1d
#define FFTWDFTR2C2D fftwf_plan_dft_r2c_2d
#define FFTWDFTC2R2D fftwf_plan_dft_c2r_2d
#define FFTWDESTROYPLAN fftwf_destroy_plan
#define FFTWIMPORTWIS fftwf_import_wisdom_from_file
#define FFTWEXECUTE fftwf_execute
#define FFTWEXECUTEDFTR2C fftwf_execute_dft_r2c

#else

#define FFTFLOAT double
#define FITSIOIMG DOUBLE_IMG
#define TFITSFLOAT TDOUBLE
#define FFTWMALLOC fftw_malloc
#define FFTWFREE fftw_free
#define FFTWCOMPLEX fftw_complex
#define FFTWPLAN fftw_plan
#define FFTWDFTR2C1D fftw_plan_dft_r2c_1d
#define FFTWDFTC2R1D fftw_plan_dft_c2r_1d
#define FFTWDFTR2C2D fftw_plan_dft_r2c_2d
#define FFTWDFTC2R2D fftw_plan_dft_c2r_2d
#define FFTWDESTROYPLAN fftw_destroy_plan
#define FFTWIMPORTWIS fftw_import_wisdom_from_file
#define FFTWEXECUTE fftw_execute
#define FFTWEXECUTEDFTR2C fftw_execute_dft_r2c

#endif

#endif
