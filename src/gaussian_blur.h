#ifndef _GAUSSIAN_BLUR_H_
#define _GAUSSIAN_BLUR_H_

#include <stdint.h>

/*******************************************************************************
 * Definitions
 ******************************************************************************/
#define FRAC_Q31(x)  ((int32_t) ((x) < 1 ? ((x) >= -1 ? (x) * 0x80000000 : 0x80000000) : 0x7FFFFFFF))
#define INV_ARR(x)   FRAC_Q31(1.0 / (x + x + 1))
  

/*******************************************************************************
 * API
 ******************************************************************************/

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus*/

void boxesForGauss(float sigma, uint32_t n, int sizes[]);

/*  
 * in_out: The buffer of input image is overwritten by the output image. 
 * tmp: The temporary buffer has same size as input buffer.
 */
void blur(uint32_t in_out[], uint32_t tmp[], uint32_t w, uint32_t h, int bxs[]);

#if defined(__cplusplus)
}
#endif /* __cplusplus*/

#endif /* _GAUSSIAN_BLUR_H_ */
