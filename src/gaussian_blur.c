#include <stdint.h>
#include <math.h>
#include "gaussian_blur.h"
#include <intrinsics.h>

/*
 * Tested with IAR EWARM
 * See http://blog.ivank.net/fastest-gaussian-blur.html
 * Author Howard Liu <howardliu7874@hotmail.com>
 */

static int32_t iarr_table[] = {
INV_ARR(0),  INV_ARR(1),  INV_ARR(2),  INV_ARR(3),  INV_ARR(4),
INV_ARR(5),  INV_ARR(6),  INV_ARR(7),  INV_ARR(8),  INV_ARR(9),
INV_ARR(10), INV_ARR(11), INV_ARR(12), INV_ARR(13), INV_ARR(14),
INV_ARR(15), INV_ARR(16), INV_ARR(17), INV_ARR(18), INV_ARR(19),
INV_ARR(20), INV_ARR(21), INV_ARR(22), INV_ARR(23), INV_ARR(24)};
 
void boxesForGauss(float sigma, uint32_t n, int sizes[])
{
    float wIdeal = sqrt((12.0*sigma*sigma/n)+1);
    int32_t wl = (int32_t) wIdeal;
    if (wl % 2 == 0) wl--;

    int32_t wu = wl + 2;

    float mIdeal = (12.0*sigma*sigma - n*wl*wl - 4*n*wl - 3*n)/(-4*wl - 4);
    int32_t m = (int32_t) (mIdeal + 0.5);

    for (int i = 0; i < n; ++i)
    {
        sizes[i] = ((i < m ? wl : wu) - 1) / 2;
    }
}

static void boxBlurH_4(uint32_t scl[], uint32_t tcl[], uint32_t w, uint32_t h, uint32_t r)
{
    int32_t iarr_q31 = iarr_table[r];
    uint32_t acc = 0x4000; // 0.5 with Q15 for round
    
    for(uint32_t row = 0; row < h; ++row)
    {
        uint32_t ti = row * w;
        uint32_t li = ti;              // left index
        uint32_t ri = ti + r;          // right index

        // first value        
        uint32_t st = scl[ti];
        uint32_t fvRB = __UXTB16(st);
        uint32_t fvAG = __UXTB16(st >> 8);

        
        // last value
        st = scl[ti + w - 1];
        uint32_t lvRB = __UXTB16(st);
        uint32_t lvAG = __UXTB16(st >> 8); 

        /* It's safe to use scaler multiplication here */
        uint32_t accAG = fvAG * (r + 1);
        uint32_t accRB = fvRB * (r + 1);

        for(uint32_t j = 0; j < r; ++j)
        {
            uint32_t s = scl[ti + j];
            
            accRB = __UXTAB16(accRB, s);
            accAG = __UXTAB16(accAG, (s >> 8));
        }

        for(uint32_t j = 0; j <= r ; ++j)
        {
            uint32_t sr = scl[ri];
            
            //acc = acc + scl[ri] - fv
            accRB = __UXTAB16(__UQSUB16(accRB, fvRB), sr);
            accAG = __UXTAB16(__UQSUB16(accAG, fvAG), (sr >> 8));
            
            //d = acc/(r+r+1)  
            uint32_t pA = __SMLAWT(iarr_q31, accAG, acc) << 1; // A component in top halfword
            uint32_t pG = __SMLAWB(iarr_q31, accAG, acc) >> 15; // G component in bottom halfword
            uint32_t pAG = __PKHTB(pA, pG, 0) << 8;
            
            uint32_t pR = __SMLAWT(iarr_q31, accRB, acc) << 1; // R component in top halfword
            uint32_t pB = __SMLAWB(iarr_q31, accRB, acc) >> 15; // B component in bottom halfword
           
            /* The high byte of each halfword in packed 32bit words
               is always zero, it's saft to OR them together. */
            tcl[ti] = __PKHTB(pR, pB, 0) | pAG; 
           
            ++ti;
            ++ri;
        }

        for(uint32_t j = r + 1; j < w - r; ++j)
        {
            //acc = acc + scl[ri] - scl[li]
            uint32_t sl = scl[li];
            uint32_t sr = scl[ri];
            
            accRB = __UXTAB16(__UQSUB16(accRB, __UXTB16(sl)), sr);
            accAG = __UXTAB16(__UQSUB16(accAG, __UXTB16(sl >> 8)), (sr >> 8));
            
            //d = acc/(r+r+1)
            uint32_t pA = __SMLAWT(iarr_q31, accAG, acc) << 1; // A component in top halfword
            uint32_t pG = __SMLAWB(iarr_q31, accAG, acc) >> 15; // G component in bottom halfword
            uint32_t pAG = __PKHTB(pA, pG, 0) << 8;
            
            uint32_t pR = __SMLAWT(iarr_q31, accRB, acc) << 1; // R component in top halfword
            uint32_t pB = __SMLAWB(iarr_q31, accRB, acc) >> 15; // B component in bottom halfword
                       
            tcl[ti] = __PKHTB(pR, pB, 0) | pAG;
            
            ++ti;
            ++ri;
            ++li;
        }

        for(uint32_t j = w - r; j < w; ++j)
        {
            uint32_t sl = scl[li];
            
            //acc = acc + lv - scl[li]
            accRB = __UQADD16(__UQSUB16(accRB, __UXTB16(sl)), lvRB);
            accAG = __UQADD16(__UQSUB16(accAG, __UXTB16(sl >> 8)), lvAG);

            //d = acc/(r+r+1)
            uint32_t pA = __SMLAWT(iarr_q31, accAG, acc) << 1; // A component in top halfword
            uint32_t pG = __SMLAWB(iarr_q31, accAG, acc) >> 15; // G component in bottom halfword
            uint32_t pAG = __PKHTB(pA, pG, 0) << 8;
            
            uint32_t pR = __SMLAWT(iarr_q31, accRB, acc) << 1; // R component in top halfword
            uint32_t pB = __SMLAWB(iarr_q31, accRB, acc) >> 15; // B component in bottom halfword
           
            tcl[ti] = __PKHTB(pR, pB, 0) | pAG;    
            
            ++ti;
            ++li;
        }
    }
}



static void boxBlurT_4(uint32_t scl[], uint32_t tcl[], uint32_t w, uint32_t h, uint32_t r)
{
    int32_t iarr_q31 = iarr_table[r];
    uint32_t acc = 0x4000; // 0.5 with Q15 for round
    
    for(uint32_t col = 0; col < w; ++col)
    {
        uint32_t ti = col;
        uint32_t li = ti;
        uint32_t ri = ti + r*w;
        
        // first value
        uint32_t st = scl[ti];
        uint32_t fvRB = __UXTB16(st);
        uint32_t fvAG = __UXTB16(st >> 8);
        
        // last value
        st = scl[ti + w*(h-1)];
        uint32_t lvRB = __UXTB16(st);
        uint32_t lvAG = __UXTB16(st >> 8); 

        /* It's safe to use scaler multiplication here */
        uint32_t accRB = fvRB * (r + 1);
        uint32_t accAG = fvAG * (r + 1);
        

        for(uint32_t j = 0; j < r; ++j)
        {
            uint32_t s = scl[ti + j * w];
            
            accRB = __UXTAB16(accRB, s);
            accAG = __UXTAB16(accAG, (s >> 8));
        }

        for(uint32_t j = 0; j <= r; ++j)
        {
            uint32_t sr = scl[ri];
            
            //acc = acc + scl[ri] - fv
            accRB = __UXTAB16(__UQSUB16(accRB, fvRB), sr);
            accAG = __UXTAB16(__UQSUB16(accAG, fvAG), (sr >> 8));

            //d = acc/(r+r+1)
            uint32_t pA = __SMLAWT(iarr_q31, accAG, acc) << 1; // A component in top halfword
            uint32_t pG = __SMLAWB(iarr_q31, accAG, acc) >> 15; // G component in bottom halfword
            uint32_t pAG = __PKHTB(pA, pG, 0) << 8;
            
            uint32_t pR = __SMLAWT(iarr_q31, accRB, acc) << 1; // R component in top halfword
            uint32_t pB = __SMLAWB(iarr_q31, accRB, acc) >> 15; // B component in bottom halfword
                       
            tcl[ti] = __PKHTB(pR, pB, 0) | pAG; 
           
            ri = ri + w;
            ti = ti + w;
        }

        for(uint32_t j = r + 1; j < h - r; ++j)
        {
            //acc = acc + scl[ri] - scl[li]
            uint32_t sl = scl[li];
            uint32_t sr = scl[ri];
            
            accRB = __UXTAB16(__UQSUB16(accRB, __UXTB16(sl)), sr);
            accAG = __UXTAB16(__UQSUB16(accAG, __UXTB16(sl >> 8)), (sr >> 8));
            
            //d = acc/(r+r+1)  
            uint32_t pA = __SMLAWT(iarr_q31, accAG, acc) << 1; // A component in top halfword
            uint32_t pG = __SMLAWB(iarr_q31, accAG, acc) >> 15; // G component in bottom halfword
            uint32_t pAG = __PKHTB(pA, pG, 0) << 8;
            
            uint32_t pR = __SMLAWT(iarr_q31, accRB, acc) << 1; // R component in top halfword
            uint32_t pB = __SMLAWB(iarr_q31, accRB, acc) >> 15; // B component in bottom halfword
                       
            tcl[ti] = __PKHTB(pR, pB, 0) | pAG;  
            
            li = li + w;
            ri = ri + w;
            ti = ti + w;
        }

        for(uint32_t j = h - r; j < h  ; ++j)
        {
            uint32_t sl = scl[li];
            
            //acc = acc + lv - scl[li]
            accRB = __UQADD16(__UQSUB16(accRB, __UXTB16(sl)), lvRB);
            accAG = __UQADD16(__UQSUB16(accAG, __UXTB16(sl >> 8)), lvAG);

            //d = acc/(r+r+1)
            uint32_t pA = __SMLAWT(iarr_q31, accAG, acc) << 1; // A component in top halfword
            uint32_t pG = __SMLAWB(iarr_q31, accAG, acc) >> 15; // G component in bottom halfword
            uint32_t pAG = __PKHTB(pA, pG, 0) << 8;
            
            uint32_t pR = __SMLAWT(iarr_q31, accRB, acc) << 1; // R component in top halfword
            uint32_t pB = __SMLAWB(iarr_q31, accRB, acc) >> 15; // B component in bottom halfword
                       
            tcl[ti] = __PKHTB(pR, pB, 0) | pAG;    
            
            li = li + w;
            ti = ti + w;
        }
    }
}

static void boxBlur_4(uint32_t in_out[], uint32_t tmp[], uint32_t w, uint32_t h, uint32_t r)
{

    boxBlurH_4(in_out, tmp, w, h, r);
    boxBlurT_4(tmp, in_out, w, h, r);
}

void blur(uint32_t in_out[], uint32_t tmp[], uint32_t w, uint32_t h, int bxs[])
{
    boxBlur_4(in_out, tmp, w, h, bxs[0]);
    boxBlur_4(in_out, tmp, w, h, bxs[1]);
    boxBlur_4(in_out, tmp, w, h, bxs[2]);
}
