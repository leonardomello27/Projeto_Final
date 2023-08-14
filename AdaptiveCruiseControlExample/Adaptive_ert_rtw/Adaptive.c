/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Adaptive.c
 *
 * Code generated for Simulink model 'Adaptive'.
 *
 * Model version                  : 11.0
 * Simulink Coder version         : 9.9 (R2023a) 19-Nov-2022
 * C/C++ source code generated on : Tue Aug  1 17:39:14 2023
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "Adaptive.h"
#include "rtwtypes.h"
#include "Adaptive_private.h"
#include <string.h>
#include <emmintrin.h>
#include <math.h>
#include "rt_nonfinite.h"

/* Named constants for MATLAB Function: '<S33>/optimizer' */
#define Adaptive_RMDscale              (0.02)
#define Adaptive_RMVscale              (0.2)
#define Adaptive_degrees               (4)
#define Adaptive_ny                    (2)
#define Adaptive_p                     (30)
#define Adaptive_voff                  (0.4)

/* Block states (default storage) */
DW_Adaptive_T Adaptive_DW;

/* External inputs (root inport signals with default storage) */
ExtU_Adaptive_T Adaptive_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_Adaptive_T Adaptive_Y;

/* Real-time model */
static RT_MODEL_Adaptive_T Adaptive_M_;
RT_MODEL_Adaptive_T *const Adaptive_M = &Adaptive_M_;

/* Forward declaration for local functions */
static real_T Adaptive_norm(const real_T x[4]);
static real_T Adaptive_maximum(const real_T x[4]);
static real_T Adaptive_xnrm2(int32_T n, const real_T x[16], int32_T ix0);
static void Adaptive_xgemv(int32_T b_m, int32_T n, const real_T b_A[16], int32_T
  ia0, const real_T x[16], int32_T ix0, real_T y[4]);
static void Adaptive_xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0,
  const real_T y[4], real_T b_A[16], int32_T ia0);
static real_T Adaptive_KWIKfactor(const real_T b_Ac[384], const int32_T iC[96],
  int32_T nA, const real_T b_Linv[16], real_T RLinv[16], real_T D[16], real_T
  b_H[16], int32_T n);
static void Adaptive_DropConstraint(int32_T kDrop, boolean_T iA[96], int32_T *nA,
  int32_T iC[96]);
static void Adaptive_qpkwik(const real_T b_Linv[16], const real_T b_Hinv[16],
  const real_T f[4], const real_T b_Ac[384], const real_T b[96], boolean_T iA[96],
  int32_T maxiter, real_T FeasTol, real_T x[4], real_T lambda[96], int32_T
  *status);

/*
 * Output and update for atomic system:
 *    '<S1>/DataTypeConversion_L0'
 *    '<S1>/DataTypeConversion_amax'
 *    '<S1>/DataTypeConversion_amin'
 *    '<S1>/DataTypeConversion_atrack'
 */
void Adaptive_DataTypeConversion_L0(real_T rtu_u, real_T *rty_y)
{
  *rty_y = rtu_u;
}

real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

/* Function for MATLAB Function: '<S33>/optimizer' */
static real_T Adaptive_norm(const real_T x[4])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = fabs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = fabs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = fabs(x[3]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * sqrt(y);
}

/* Function for MATLAB Function: '<S33>/optimizer' */
static real_T Adaptive_maximum(const real_T x[4])
{
  real_T ex;
  int32_T idx;
  int32_T k;
  if (!rtIsNaN(x[0])) {
    idx = 1;
  } else {
    boolean_T exitg1;
    idx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 5)) {
      if (!rtIsNaN(x[k - 1])) {
        idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (idx == 0) {
    ex = x[0];
  } else {
    ex = x[idx - 1];
    for (k = idx + 1; k < 5; k++) {
      real_T x_0;
      x_0 = x[k - 1];
      if (ex < x_0) {
        ex = x_0;
      }
    }
  }

  return ex;
}

/* Function for MATLAB Function: '<S33>/optimizer' */
static real_T Adaptive_xnrm2(int32_T n, const real_T x[16], int32_T ix0)
{
  real_T y;
  int32_T k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        real_T absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T a;
  real_T b;
  real_T y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = sqrt(a * a + 1.0) * b;
  } else if (a > b) {
    b /= a;
    y = sqrt(b * b + 1.0) * a;
  } else if (rtIsNaN(b)) {
    y = (rtNaN);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

/* Function for MATLAB Function: '<S33>/optimizer' */
static void Adaptive_xgemv(int32_T b_m, int32_T n, const real_T b_A[16], int32_T
  ia0, const real_T x[16], int32_T ix0, real_T y[4])
{
  int32_T b_iy;
  int32_T ia;
  if ((b_m != 0) && (n != 0)) {
    int32_T b;
    if (n - 1 >= 0) {
      memset(&y[0], 0, (uint32_T)n * sizeof(real_T));
    }

    b = ((n - 1) << 2) + ia0;
    for (b_iy = ia0; b_iy <= b; b_iy += 4) {
      real_T c;
      int32_T d;
      c = 0.0;
      d = b_iy + b_m;
      for (ia = b_iy; ia < d; ia++) {
        c += x[((ix0 + ia) - b_iy) - 1] * b_A[ia - 1];
      }

      ia = (b_iy - ia0) >> 2;
      y[ia] += c;
    }
  }
}

/* Function for MATLAB Function: '<S33>/optimizer' */
static void Adaptive_xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0,
  const real_T y[4], real_T b_A[16], int32_T ia0)
{
  int32_T ijA;
  int32_T j;
  if (!(alpha1 == 0.0)) {
    int32_T jA;
    jA = ia0;
    for (j = 0; j < n; j++) {
      real_T temp;
      temp = y[j];
      if (temp != 0.0) {
        int32_T b;
        temp *= alpha1;
        b = b_m + jA;
        for (ijA = jA; ijA < b; ijA++) {
          b_A[ijA - 1] += b_A[((ix0 + ijA) - jA) - 1] * temp;
        }
      }

      jA += 4;
    }
  }
}

/* Function for MATLAB Function: '<S33>/optimizer' */
static real_T Adaptive_KWIKfactor(const real_T b_Ac[384], const int32_T iC[96],
  int32_T nA, const real_T b_Linv[16], real_T RLinv[16], real_T D[16], real_T
  b_H[16], int32_T n)
{
  __m128d tmp;
  real_T Q[16];
  real_T R[16];
  real_T TL[16];
  real_T b_A[16];
  real_T tau[4];
  real_T work[4];
  real_T Status;
  real_T atmp;
  real_T b_A_0;
  real_T beta1;
  int32_T b_lastv;
  int32_T c_lastc;
  int32_T exitg1;
  int32_T g;
  int32_T h_k;
  int32_T ii;
  int32_T k_i;
  int32_T knt;
  int32_T vectorUB_tmp;
  boolean_T exitg2;
  Status = 1.0;
  memset(&RLinv[0], 0, sizeof(real_T) << 4U);
  for (ii = 0; ii < nA; ii++) {
    b_lastv = iC[ii];
    for (k_i = 0; k_i <= 2; k_i += 2) {
      _mm_storeu_pd(&RLinv[k_i + (ii << 2)], _mm_add_pd(_mm_add_pd(_mm_add_pd
        (_mm_mul_pd(_mm_set1_pd(b_Ac[b_lastv - 1]), _mm_loadu_pd(&b_Linv[k_i])),
         _mm_mul_pd(_mm_loadu_pd(&b_Linv[k_i + 4]), _mm_set1_pd(b_Ac[b_lastv +
        95]))), _mm_mul_pd(_mm_loadu_pd(&b_Linv[k_i + 8]), _mm_set1_pd
                           (b_Ac[b_lastv + 191]))), _mm_mul_pd(_mm_loadu_pd
        (&b_Linv[k_i + 12]), _mm_set1_pd(b_Ac[b_lastv + 287]))));
    }
  }

  memcpy(&b_A[0], &RLinv[0], sizeof(real_T) << 4U);
  tau[0] = 0.0;
  work[0] = 0.0;
  tau[1] = 0.0;
  work[1] = 0.0;
  tau[2] = 0.0;
  work[2] = 0.0;
  tau[3] = 0.0;
  work[3] = 0.0;
  for (k_i = 0; k_i < 4; k_i++) {
    ii = (k_i << 2) + k_i;
    if (k_i + 1 < 4) {
      atmp = b_A[ii];
      b_lastv = ii + 2;
      tau[k_i] = 0.0;
      beta1 = Adaptive_xnrm2(3 - k_i, b_A, ii + 2);
      if (beta1 != 0.0) {
        b_A_0 = b_A[ii];
        beta1 = rt_hypotd_snf(b_A_0, beta1);
        if (b_A_0 >= 0.0) {
          beta1 = -beta1;
        }

        if (fabs(beta1) < 1.0020841800044864E-292) {
          knt = 0;
          h_k = (ii - k_i) + 4;
          do {
            knt++;
            g = (((((h_k - ii) - 1) / 2) << 1) + ii) + 2;
            vectorUB_tmp = g - 2;
            for (c_lastc = b_lastv; c_lastc <= vectorUB_tmp; c_lastc += 2) {
              tmp = _mm_loadu_pd(&b_A[c_lastc - 1]);
              _mm_storeu_pd(&b_A[c_lastc - 1], _mm_mul_pd(tmp, _mm_set1_pd
                (9.9792015476736E+291)));
            }

            for (c_lastc = g; c_lastc <= h_k; c_lastc++) {
              b_A[c_lastc - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while ((fabs(beta1) < 1.0020841800044864E-292) && (knt < 20));

          beta1 = rt_hypotd_snf(atmp, Adaptive_xnrm2(3 - k_i, b_A, ii + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[k_i] = (beta1 - atmp) / beta1;
          atmp = 1.0 / (atmp - beta1);
          for (c_lastc = b_lastv; c_lastc <= vectorUB_tmp; c_lastc += 2) {
            tmp = _mm_loadu_pd(&b_A[c_lastc - 1]);
            _mm_storeu_pd(&b_A[c_lastc - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
          }

          for (c_lastc = g; c_lastc <= h_k; c_lastc++) {
            b_A[c_lastc - 1] *= atmp;
          }

          for (b_lastv = 0; b_lastv < knt; b_lastv++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[k_i] = (beta1 - b_A_0) / beta1;
          atmp = 1.0 / (b_A_0 - beta1);
          g = (ii - k_i) + 4;
          c_lastc = (((((g - ii) - 1) / 2) << 1) + ii) + 2;
          knt = c_lastc - 2;
          for (h_k = b_lastv; h_k <= knt; h_k += 2) {
            tmp = _mm_loadu_pd(&b_A[h_k - 1]);
            _mm_storeu_pd(&b_A[h_k - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
          }

          for (h_k = c_lastc; h_k <= g; h_k++) {
            b_A[h_k - 1] *= atmp;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = 1.0;
      if (tau[k_i] != 0.0) {
        b_lastv = 4 - k_i;
        c_lastc = (ii - k_i) + 3;
        while ((b_lastv > 0) && (b_A[c_lastc] == 0.0)) {
          b_lastv--;
          c_lastc--;
        }

        c_lastc = 3 - k_i;
        exitg2 = false;
        while ((!exitg2) && (c_lastc > 0)) {
          knt = (((c_lastc - 1) << 2) + ii) + 4;
          h_k = knt;
          do {
            exitg1 = 0;
            if (h_k + 1 <= knt + b_lastv) {
              if (b_A[h_k] != 0.0) {
                exitg1 = 1;
              } else {
                h_k++;
              }
            } else {
              c_lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        b_lastv = 0;
        c_lastc = 0;
      }

      if (b_lastv > 0) {
        Adaptive_xgemv(b_lastv, c_lastc, b_A, ii + 5, b_A, ii + 1, work);
        Adaptive_xgerc(b_lastv, c_lastc, -tau[k_i], ii + 1, work, b_A, ii + 5);
      }

      b_A[ii] = atmp;
    } else {
      tau[3] = 0.0;
    }
  }

  for (k_i = 0; k_i < 4; k_i++) {
    for (ii = 0; ii <= k_i; ii++) {
      b_lastv = k_i << 2;
      R[ii + b_lastv] = b_A[b_lastv + ii];
    }

    for (ii = k_i + 2; ii < 5; ii++) {
      R[(ii + (k_i << 2)) - 1] = 0.0;
    }

    work[k_i] = 0.0;
  }

  for (k_i = 3; k_i >= 0; k_i--) {
    ii = ((k_i << 2) + k_i) + 5;
    if (k_i + 1 < 4) {
      b_A[ii - 5] = 1.0;
      if (tau[k_i] != 0.0) {
        b_lastv = 4 - k_i;
        c_lastc = ii - k_i;
        while ((b_lastv > 0) && (b_A[c_lastc - 2] == 0.0)) {
          b_lastv--;
          c_lastc--;
        }

        c_lastc = 3 - k_i;
        exitg2 = false;
        while ((!exitg2) && (c_lastc > 0)) {
          knt = ((c_lastc - 1) << 2) + ii;
          h_k = knt;
          do {
            exitg1 = 0;
            if (h_k <= (knt + b_lastv) - 1) {
              if (b_A[h_k - 1] != 0.0) {
                exitg1 = 1;
              } else {
                h_k++;
              }
            } else {
              c_lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        b_lastv = 0;
        c_lastc = 0;
      }

      if (b_lastv > 0) {
        Adaptive_xgemv(b_lastv, c_lastc, b_A, ii, b_A, ii - 4, work);
        Adaptive_xgerc(b_lastv, c_lastc, -tau[k_i], ii - 4, work, b_A, ii);
      }

      h_k = (ii - k_i) - 1;
      c_lastc = (((((h_k - ii) + 4) / 2) << 1) + ii) - 3;
      knt = c_lastc - 2;
      for (b_lastv = ii - 3; b_lastv <= knt; b_lastv += 2) {
        tmp = _mm_loadu_pd(&b_A[b_lastv - 1]);
        _mm_storeu_pd(&b_A[b_lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(-tau[k_i])));
      }

      for (b_lastv = c_lastc; b_lastv <= h_k; b_lastv++) {
        b_A[b_lastv - 1] *= -tau[k_i];
      }
    }

    b_A[ii - 5] = 1.0 - tau[k_i];
    for (b_lastv = 0; b_lastv < k_i; b_lastv++) {
      b_A[(ii - b_lastv) - 6] = 0.0;
    }
  }

  for (k_i = 0; k_i < 4; k_i++) {
    ii = k_i << 2;
    Q[ii] = b_A[ii];
    Q[ii + 1] = b_A[ii + 1];
    Q[ii + 2] = b_A[ii + 2];
    Q[ii + 3] = b_A[ii + 3];
  }

  k_i = 0;
  do {
    exitg1 = 0;
    if (k_i <= nA - 1) {
      if (fabs(R[(k_i << 2) + k_i]) < 1.0E-12) {
        Status = -2.0;
        exitg1 = 1;
      } else {
        k_i++;
      }
    } else {
      for (k_i = 0; k_i < n; k_i++) {
        for (ii = 0; ii < n; ii++) {
          b_lastv = k_i << 2;
          c_lastc = ii << 2;
          TL[k_i + c_lastc] = ((b_Linv[b_lastv + 1] * Q[c_lastc + 1] +
                                b_Linv[b_lastv] * Q[c_lastc]) + b_Linv[b_lastv +
                               2] * Q[c_lastc + 2]) + b_Linv[b_lastv + 3] *
            Q[c_lastc + 3];
        }
      }

      memset(&RLinv[0], 0, sizeof(real_T) << 4U);
      for (b_lastv = nA; b_lastv >= 1; b_lastv--) {
        k_i = (b_lastv - 1) << 2;
        ii = (b_lastv + k_i) - 1;
        RLinv[ii] = 1.0;
        for (c_lastc = b_lastv; c_lastc <= nA; c_lastc++) {
          h_k = (((c_lastc - 1) << 2) + b_lastv) - 1;
          RLinv[h_k] /= R[ii];
        }

        if (b_lastv > 1) {
          for (c_lastc = 0; c_lastc <= b_lastv - 2; c_lastc++) {
            for (knt = b_lastv; knt <= nA; knt++) {
              ii = (knt - 1) << 2;
              h_k = ii + c_lastc;
              RLinv[h_k] -= RLinv[(ii + b_lastv) - 1] * R[k_i + c_lastc];
            }
          }
        }
      }

      for (b_lastv = 0; b_lastv < n; b_lastv++) {
        for (c_lastc = b_lastv + 1; c_lastc <= n; c_lastc++) {
          k_i = ((c_lastc - 1) << 2) + b_lastv;
          b_H[k_i] = 0.0;
          for (knt = nA + 1; knt <= n; knt++) {
            ii = (knt - 1) << 2;
            b_H[k_i] -= TL[(ii + c_lastc) - 1] * TL[ii + b_lastv];
          }

          b_H[(c_lastc + (b_lastv << 2)) - 1] = b_H[k_i];
        }
      }

      for (b_lastv = 0; b_lastv < nA; b_lastv++) {
        for (c_lastc = 0; c_lastc < n; c_lastc++) {
          k_i = (b_lastv << 2) + c_lastc;
          D[k_i] = 0.0;
          for (knt = b_lastv + 1; knt <= nA; knt++) {
            ii = (knt - 1) << 2;
            D[k_i] += TL[ii + c_lastc] * RLinv[ii + b_lastv];
          }
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return Status;
}

/* Function for MATLAB Function: '<S33>/optimizer' */
static void Adaptive_DropConstraint(int32_T kDrop, boolean_T iA[96], int32_T *nA,
  int32_T iC[96])
{
  int32_T i;
  if (kDrop > 0) {
    iA[iC[kDrop - 1] - 1] = false;
    if (kDrop < *nA) {
      int32_T b;
      b = *nA - 1;
      for (i = kDrop; i <= b; i++) {
        iC[i - 1] = iC[i];
      }
    }

    iC[*nA - 1] = 0;
    (*nA)--;
  }
}

/* Function for MATLAB Function: '<S33>/optimizer' */
static void Adaptive_qpkwik(const real_T b_Linv[16], const real_T b_Hinv[16],
  const real_T f[4], const real_T b_Ac[384], const real_T b[96], boolean_T iA[96],
  int32_T maxiter, real_T FeasTol, real_T x[4], real_T lambda[96], int32_T
  *status)
{
  __m128d tmp_0;
  __m128d tmp_1;
  __m128d tmp_2;
  __m128d tmp_3;
  real_T cTol[96];
  real_T D[16];
  real_T RLinv[16];
  real_T U[16];
  real_T b_H[16];
  real_T Opt[8];
  real_T Rhs[8];
  real_T r[4];
  real_T z[4];
  real_T Xnorm0;
  real_T cMin;
  real_T cVal;
  real_T rMin;
  real_T t;
  real_T zTa;
  real_T zTa_tmp;
  int32_T iC[96];
  int32_T U_tmp;
  int32_T U_tmp_0;
  int32_T b_exponent;
  int32_T exitg1;
  int32_T exitg3;
  int32_T exponent;
  int32_T i;
  int32_T iSave;
  int32_T nA;
  int32_T r_tmp;
  int32_T tmp;
  boolean_T ColdReset;
  boolean_T DualFeasible;
  boolean_T cTolComputed;
  boolean_T exitg2;
  boolean_T exitg4;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = 0.0;
  *status = 1;
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  r[3] = 0.0;
  rMin = 0.0;
  cTolComputed = false;
  for (i = 0; i < 96; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (i = 0; i < 96; i++) {
    if (iA[i]) {
      nA++;
      iC[nA - 1] = i + 1;
    }
  }

  guard1 = false;
  if (nA > 0) {
    memset(&Opt[0], 0, sizeof(real_T) << 3U);
    Rhs[0] = f[0];
    Rhs[4] = 0.0;
    Rhs[1] = f[1];
    Rhs[5] = 0.0;
    Rhs[2] = f[2];
    Rhs[6] = 0.0;
    Rhs[3] = f[3];
    Rhs[7] = 0.0;
    DualFeasible = false;
    tmp = (int32_T)rt_roundd_snf(0.3 * (real_T)nA);
    ColdReset = false;
    do {
      exitg3 = 0;
      if ((!DualFeasible) && (nA > 0) && (*status <= maxiter)) {
        Xnorm0 = Adaptive_KWIKfactor(b_Ac, iC, nA, b_Linv, RLinv, D, b_H,
          Adaptive_degrees);
        if (Xnorm0 < 0.0) {
          if (ColdReset) {
            *status = -2;
            exitg3 = 2;
          } else {
            nA = 0;
            memset(&iA[0], 0, 96U * sizeof(boolean_T));
            memset(&iC[0], 0, 96U * sizeof(int32_T));
            ColdReset = true;
          }
        } else {
          for (i = 0; i < nA; i++) {
            Rhs[i + 4] = b[iC[i] - 1];
            for (r_tmp = i + 1; r_tmp <= nA; r_tmp++) {
              U_tmp_0 = ((i << 2) + r_tmp) - 1;
              U[U_tmp_0] = 0.0;
              for (iSave = 0; iSave < nA; iSave++) {
                U_tmp = iSave << 2;
                U[U_tmp_0] += RLinv[(U_tmp + r_tmp) - 1] * RLinv[U_tmp + i];
              }

              U[i + ((r_tmp - 1) << 2)] = U[U_tmp_0];
            }
          }

          for (i = 0; i < 4; i++) {
            Opt[i] = ((b_H[i + 4] * Rhs[1] + b_H[i] * Rhs[0]) + b_H[i + 8] *
                      Rhs[2]) + b_H[i + 12] * Rhs[3];
            for (r_tmp = 0; r_tmp < nA; r_tmp++) {
              Opt[i] += D[(r_tmp << 2) + i] * Rhs[r_tmp + 4];
            }
          }

          Xnorm0 = -1.0E-12;
          i = -1;
          for (r_tmp = 0; r_tmp < nA; r_tmp++) {
            iSave = r_tmp << 2;
            Opt[r_tmp + 4] = ((D[iSave + 1] * Rhs[1] + D[iSave] * Rhs[0]) +
                              D[iSave + 2] * Rhs[2]) + D[iSave + 3] * Rhs[3];
            for (iSave = 0; iSave < nA; iSave++) {
              Opt[r_tmp + 4] += U[(iSave << 2) + r_tmp] * Rhs[iSave + 4];
            }

            cMin = Opt[r_tmp + 4];
            lambda[iC[r_tmp] - 1] = cMin;
            if ((cMin < Xnorm0) && (r_tmp + 1 <= nA)) {
              i = r_tmp;
              Xnorm0 = cMin;
            }
          }

          if (i + 1 <= 0) {
            DualFeasible = true;
            x[0] = Opt[0];
            x[1] = Opt[1];
            x[2] = Opt[2];
            x[3] = Opt[3];
          } else {
            (*status)++;
            if (tmp <= 5) {
              r_tmp = 5;
            } else {
              r_tmp = tmp;
            }

            if (*status > r_tmp) {
              nA = 0;
              memset(&iA[0], 0, 96U * sizeof(boolean_T));
              memset(&iC[0], 0, 96U * sizeof(int32_T));
              ColdReset = true;
            } else {
              lambda[iC[i] - 1] = 0.0;
              Adaptive_DropConstraint(i + 1, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          memset(&lambda[0], 0, 96U * sizeof(real_T));
          Xnorm0 = f[1];
          cMin = f[0];
          cVal = f[2];
          t = f[3];
          for (tmp = 0; tmp <= 2; tmp += 2) {
            tmp_3 = _mm_set1_pd(-1.0);
            _mm_storeu_pd(&x[tmp], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
              (_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 4]), tmp_3), _mm_set1_pd
               (Xnorm0)), _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp]),
              tmp_3), _mm_set1_pd(cMin))), _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd
              (&b_Hinv[tmp + 8]), tmp_3), _mm_set1_pd(cVal))), _mm_mul_pd
              (_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 12]), tmp_3), _mm_set1_pd(t))));
          }
        }

        exitg3 = 1;
      }
    } while (exitg3 == 0);

    if (exitg3 == 1) {
      guard1 = true;
    }
  } else {
    Xnorm0 = f[1];
    cMin = f[0];
    cVal = f[2];
    t = f[3];
    for (tmp = 0; tmp <= 2; tmp += 2) {
      tmp_3 = _mm_set1_pd(-1.0);
      _mm_storeu_pd(&x[tmp], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
        (_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 4]), tmp_3), _mm_set1_pd(Xnorm0)),
        _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp]), tmp_3), _mm_set1_pd
                   (cMin))), _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 8]),
        tmp_3), _mm_set1_pd(cVal))), _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd
        (&b_Hinv[tmp + 12]), tmp_3), _mm_set1_pd(t))));
    }

    guard1 = true;
  }

  if (guard1) {
    Xnorm0 = Adaptive_norm(x);
    exitg2 = false;
    while ((!exitg2) && (*status <= maxiter)) {
      cMin = -FeasTol;
      tmp = -1;
      for (i = 0; i < 96; i++) {
        if (!cTolComputed) {
          z[0] = fabs(b_Ac[i] * x[0]);
          z[1] = fabs(b_Ac[i + 96] * x[1]);
          z[2] = fabs(b_Ac[i + 192] * x[2]);
          z[3] = fabs(b_Ac[i + 288] * x[3]);
          cTol[i] = fmax(cTol[i], Adaptive_maximum(z));
        }

        if (!iA[i]) {
          cVal = ((((b_Ac[i + 96] * x[1] + b_Ac[i] * x[0]) + b_Ac[i + 192] * x[2])
                   + b_Ac[i + 288] * x[3]) - b[i]) / cTol[i];
          if (cVal < cMin) {
            cMin = cVal;
            tmp = i;
          }
        }
      }

      cTolComputed = true;
      if (tmp + 1 <= 0) {
        exitg2 = true;
      } else if (*status == maxiter) {
        *status = 0;
        exitg2 = true;
      } else {
        do {
          exitg1 = 0;
          if ((tmp + 1 > 0) && (*status <= maxiter)) {
            guard2 = false;
            if (nA == 0) {
              for (r_tmp = 0; r_tmp <= 2; r_tmp += 2) {
                _mm_storeu_pd(&z[r_tmp], _mm_add_pd(_mm_add_pd(_mm_add_pd
                  (_mm_mul_pd(_mm_loadu_pd(&b_Hinv[r_tmp + 4]), _mm_set1_pd
                              (b_Ac[tmp + 96])), _mm_mul_pd(_mm_loadu_pd
                  (&b_Hinv[r_tmp]), _mm_set1_pd(b_Ac[tmp]))), _mm_mul_pd
                  (_mm_loadu_pd(&b_Hinv[r_tmp + 8]), _mm_set1_pd(b_Ac[tmp + 192]))),
                  _mm_mul_pd(_mm_loadu_pd(&b_Hinv[r_tmp + 12]), _mm_set1_pd
                             (b_Ac[tmp + 288]))));
              }

              guard2 = true;
            } else {
              cMin = Adaptive_KWIKfactor(b_Ac, iC, nA, b_Linv, RLinv, D, b_H,
                Adaptive_degrees);
              if (cMin <= 0.0) {
                *status = -2;
                exitg1 = 1;
              } else {
                for (r_tmp = 0; r_tmp <= 14; r_tmp += 2) {
                  tmp_3 = _mm_loadu_pd(&b_H[r_tmp]);
                  _mm_storeu_pd(&U[r_tmp], _mm_mul_pd(tmp_3, _mm_set1_pd(-1.0)));
                }

                for (r_tmp = 0; r_tmp <= 2; r_tmp += 2) {
                  tmp_3 = _mm_loadu_pd(&U[r_tmp + 4]);
                  tmp_0 = _mm_loadu_pd(&U[r_tmp]);
                  tmp_1 = _mm_loadu_pd(&U[r_tmp + 8]);
                  tmp_2 = _mm_loadu_pd(&U[r_tmp + 12]);
                  _mm_storeu_pd(&z[r_tmp], _mm_add_pd(_mm_add_pd(_mm_add_pd
                    (_mm_mul_pd(tmp_3, _mm_set1_pd(b_Ac[tmp + 96])), _mm_mul_pd
                     (tmp_0, _mm_set1_pd(b_Ac[tmp]))), _mm_mul_pd(tmp_1,
                    _mm_set1_pd(b_Ac[tmp + 192]))), _mm_mul_pd(tmp_2,
                    _mm_set1_pd(b_Ac[tmp + 288]))));
                }

                for (i = 0; i < nA; i++) {
                  r_tmp = i << 2;
                  r[i] = ((D[r_tmp + 1] * b_Ac[tmp + 96] + D[r_tmp] * b_Ac[tmp])
                          + D[r_tmp + 2] * b_Ac[tmp + 192]) + D[r_tmp + 3] *
                    b_Ac[tmp + 288];
                }

                guard2 = true;
              }
            }

            if (guard2) {
              i = 0;
              cMin = 0.0;
              DualFeasible = true;
              ColdReset = true;
              if (nA > 0) {
                r_tmp = 0;
                exitg4 = false;
                while ((!exitg4) && (r_tmp <= nA - 1)) {
                  if (r[r_tmp] >= 1.0E-12) {
                    ColdReset = false;
                    exitg4 = true;
                  } else {
                    r_tmp++;
                  }
                }
              }

              if ((nA != 0) && (!ColdReset)) {
                for (r_tmp = 0; r_tmp < nA; r_tmp++) {
                  cVal = r[r_tmp];
                  if (cVal > 1.0E-12) {
                    cVal = lambda[iC[r_tmp] - 1] / cVal;
                    if ((i == 0) || (cVal < rMin)) {
                      rMin = cVal;
                      i = r_tmp + 1;
                    }
                  }
                }

                if (i > 0) {
                  cMin = rMin;
                  DualFeasible = false;
                }
              }

              cVal = b_Ac[tmp + 96];
              t = b_Ac[tmp + 192];
              zTa_tmp = b_Ac[tmp + 288];
              zTa = ((cVal * z[1] + z[0] * b_Ac[tmp]) + t * z[2]) + zTa_tmp * z
                [3];
              if (zTa <= 0.0) {
                cVal = 0.0;
                ColdReset = true;
              } else {
                cVal = (b[tmp] - (((cVal * x[1] + b_Ac[tmp] * x[0]) + t * x[2])
                                  + zTa_tmp * x[3])) / zTa;
                ColdReset = false;
              }

              if (DualFeasible && ColdReset) {
                *status = -1;
                exitg1 = 1;
              } else {
                if (ColdReset) {
                  t = cMin;
                } else if (DualFeasible) {
                  t = cVal;
                } else if (cMin < cVal) {
                  t = cMin;
                } else {
                  t = cVal;
                }

                for (r_tmp = 0; r_tmp < nA; r_tmp++) {
                  iSave = iC[r_tmp];
                  lambda[iSave - 1] -= t * r[r_tmp];
                  if ((iSave <= 96) && (lambda[iSave - 1] < 0.0)) {
                    lambda[iSave - 1] = 0.0;
                  }
                }

                lambda[tmp] += t;
                frexp(1.0, &exponent);
                if (fabs(t - cMin) < 2.2204460492503131E-16) {
                  Adaptive_DropConstraint(i, iA, &nA, iC);
                }

                if (!ColdReset) {
                  x[0] += t * z[0];
                  x[1] += t * z[1];
                  x[2] += t * z[2];
                  x[3] += t * z[3];
                  frexp(1.0, &b_exponent);
                  if (fabs(t - cVal) < 2.2204460492503131E-16) {
                    if (nA == Adaptive_degrees) {
                      *status = -1;
                      exitg1 = 1;
                    } else {
                      nA++;
                      iC[nA - 1] = tmp + 1;
                      i = nA - 1;
                      exitg4 = false;
                      while ((!exitg4) && (i + 1 > 1)) {
                        r_tmp = iC[i - 1];
                        if (iC[i] > r_tmp) {
                          exitg4 = true;
                        } else {
                          iSave = iC[i];
                          iC[i] = r_tmp;
                          iC[i - 1] = iSave;
                          i--;
                        }
                      }

                      iA[tmp] = true;
                      tmp = -1;
                      (*status)++;
                    }
                  } else {
                    (*status)++;
                  }
                } else {
                  (*status)++;
                }
              }
            }
          } else {
            cMin = Adaptive_norm(x);
            if (fabs(cMin - Xnorm0) > 0.001) {
              Xnorm0 = cMin;
              for (tmp = 0; tmp < 96; tmp++) {
                cTol[tmp] = fmax(fabs(b[tmp]), 1.0);
              }

              cTolComputed = false;
            }

            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }
  }
}

/* Model step function */
void Adaptive_step(void)
{
  real_T Bc[96];
  real_T a__1[96];
  real_T vseq[62];
  real_T rseq[60];
  real_T f[4];
  real_T xk[4];
  real_T y_innov[2];
  real_T ymax_incr[2];
  real_T ymin_incr[2];
  boolean_T ymax_incr_flag[2];
  boolean_T ymin_incr_flag[2];
  boolean_T b_Del_Save_Flag0;
  boolean_T umax_incr_flag;
  boolean_T umin_incr_flag;
  static const real_T a[8] = { -0.28835920047412961, 0.59465332667199566,
    3.6967671444929451, 0.00894917948942857, 0.00779341054829927,
    -0.056584183908856341, -0.016865406539239628, 0.094687375766569562 };

  static const real_T b_a[8] = { -4.4692446269479357E-18, -0.0013949226148762242,
    3.5355336307537657, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T c_a[4] = { -1.1168314590257729, -0.63245504018285259,
    -0.0013713976378484983, 0.0 };

  static const real_T d_a[16] = { 0.81873075307798171, 7.9960268863730558E-6,
    -0.020266516340000176, 0.0, -0.090634616406701243, 0.99998194235173521,
    0.04576843335426678, 0.0, -3.5759319333468191E-5, -7.1245318434864537E-9,
    1.0000180576482651, 0.0, 0.0, 0.0, 0.0, 1.0 };

  static const real_T e_a[8] = { -0.28621232262067559, 0.59466637961177371,
    3.6636835117083622, 0.0089491794894285132, 0.0032542098114171168,
    -0.056585231827053016, -0.014209381043202364, 0.094687375766569618 };

  static const real_T f_a[8] = { -2.7438077476480415E-20, 0.0063245553203367215,
    -1.1159380919009811E-5, -0.012649109656163445, 0.028284269046030115,
    -4.99062686434405E-6, 0.0, 1.0 };

  __m128d tmp;
  real_T rtb_xest[4];
  real_T rtb_TmpSignalConversionAtSFun_e[2];
  real_T rtb_y;
  real_T rtb_y_b;
  real_T rtb_y_d;
  real_T rtb_y_ds;
  real_T rtb_y_g;
  real_T rtb_y_j;
  real_T rtb_y_ov;
  real_T xk_0;
  real_T y_innov_0;
  real_T y_innov_1;
  int32_T d_i;
  int32_T i;
  int8_T rtb_TmpSignalConversionAtSFu_mm[2];
  uint8_T b_Mrows;
  static const real_T b_Mlim[96] = { 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
    0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
    0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4,
    0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4,
    0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4,
    0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4,
    0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.56, 0.4, 0.4, 0.4, 0.4, 0.6,
    0.6, 0.6 };

  static const real_T b_Mv[5952] = { 1.6940658945086007E-20,
    2.0328790734103208E-20, 2.0328790734103208E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 3.0493186101154812E-20,
    2.7105054312137611E-20, 3.3881317890172014E-20, 3.3881317890172014E-20,
    3.0493186101154812E-20, 3.7269449679189215E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    3.7269449679189215E-20, 4.4045713257223618E-20, 4.7433845046240819E-20,
    4.7433845046240819E-20, 5.082197683525802E-20, 5.7598240413292423E-20,
    5.7598240413292423E-20, 6.4374503991326826E-20, 6.4374503991326826E-20,
    6.4374503991326826E-20, 6.7762635780344027E-20, 6.7762635780344027E-20,
    7.1150767569361228E-20, 7.453889935837843E-20, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000007, -3.0493186101154812E-20,
    0.1000000000000001, -2.7105054312137611E-20, 0.1000000000000001,
    -3.3881317890172014E-20, 0.10000000000000012, -3.3881317890172014E-20,
    0.10000000000000012, -3.0493186101154812E-20, 0.10000000000000013,
    -3.7269449679189215E-20, 0.10000000000000013, -4.0657581468206416E-20,
    0.10000000000000014, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000016, -4.0657581468206416E-20,
    0.10000000000000017, -3.7269449679189215E-20, 0.10000000000000019,
    -4.4045713257223618E-20, 0.10000000000000019, -4.7433845046240819E-20,
    0.10000000000000021, -4.7433845046240819E-20, 0.10000000000000021,
    -5.082197683525802E-20, 0.10000000000000023, -5.7598240413292423E-20,
    0.10000000000000024, -5.7598240413292423E-20, 0.10000000000000026,
    -6.4374503991326826E-20, 0.10000000000000026, -6.4374503991326826E-20,
    0.10000000000000026, -6.4374503991326826E-20, 0.10000000000000028,
    -6.7762635780344027E-20, 0.10000000000000028, -6.7762635780344027E-20,
    0.1000000000000003, -7.1150767569361228E-20, 0.1000000000000003,
    -7.453889935837843E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, 1.6940658945086007E-20,
    2.0328790734103208E-20, 2.0328790734103208E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 3.0493186101154812E-20,
    2.7105054312137611E-20, 3.3881317890172014E-20, 3.3881317890172014E-20,
    3.0493186101154812E-20, 3.7269449679189215E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    3.7269449679189215E-20, 4.4045713257223618E-20, 4.7433845046240819E-20,
    4.7433845046240819E-20, 5.082197683525802E-20, 5.7598240413292423E-20,
    5.7598240413292423E-20, 6.4374503991326826E-20, 6.4374503991326826E-20,
    6.4374503991326826E-20, 6.7762635780344027E-20, 6.7762635780344027E-20,
    7.1150767569361228E-20, 0.0, 0.0, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000007, -3.0493186101154812E-20,
    0.1000000000000001, -2.7105054312137611E-20, 0.1000000000000001,
    -3.3881317890172014E-20, 0.10000000000000012, -3.3881317890172014E-20,
    0.10000000000000012, -3.0493186101154812E-20, 0.10000000000000013,
    -3.7269449679189215E-20, 0.10000000000000013, -4.0657581468206416E-20,
    0.10000000000000014, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000016, -4.0657581468206416E-20,
    0.10000000000000017, -3.7269449679189215E-20, 0.10000000000000019,
    -4.4045713257223618E-20, 0.10000000000000019, -4.7433845046240819E-20,
    0.10000000000000021, -4.7433845046240819E-20, 0.10000000000000021,
    -5.082197683525802E-20, 0.10000000000000023, -5.7598240413292423E-20,
    0.10000000000000024, -5.7598240413292423E-20, 0.10000000000000026,
    -6.4374503991326826E-20, 0.10000000000000026, -6.4374503991326826E-20,
    0.10000000000000026, -6.4374503991326826E-20, 0.10000000000000028,
    -6.7762635780344027E-20, 0.10000000000000028, -6.7762635780344027E-20,
    0.1000000000000003, -7.1150767569361228E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0,
    1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 3.3881317890172014E-20,
    3.3881317890172014E-20, 3.0493186101154812E-20, 3.7269449679189215E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 3.7269449679189215E-20, 4.4045713257223618E-20,
    4.7433845046240819E-20, 4.7433845046240819E-20, 5.082197683525802E-20,
    5.7598240413292423E-20, 5.7598240413292423E-20, 6.4374503991326826E-20,
    6.4374503991326826E-20, 6.4374503991326826E-20, 6.7762635780344027E-20,
    6.7762635780344027E-20, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000007, -3.0493186101154812E-20,
    0.1000000000000001, -2.7105054312137611E-20, 0.1000000000000001,
    -3.3881317890172014E-20, 0.10000000000000012, -3.3881317890172014E-20,
    0.10000000000000012, -3.0493186101154812E-20, 0.10000000000000013,
    -3.7269449679189215E-20, 0.10000000000000013, -4.0657581468206416E-20,
    0.10000000000000014, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000016, -4.0657581468206416E-20,
    0.10000000000000017, -3.7269449679189215E-20, 0.10000000000000019,
    -4.4045713257223618E-20, 0.10000000000000019, -4.7433845046240819E-20,
    0.10000000000000021, -4.7433845046240819E-20, 0.10000000000000021,
    -5.082197683525802E-20, 0.10000000000000023, -5.7598240413292423E-20,
    0.10000000000000024, -5.7598240413292423E-20, 0.10000000000000026,
    -6.4374503991326826E-20, 0.10000000000000026, -6.4374503991326826E-20,
    0.10000000000000026, -6.4374503991326826E-20, 0.10000000000000028,
    -6.7762635780344027E-20, 0.10000000000000028, -6.7762635780344027E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20,
    2.0328790734103208E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 3.0493186101154812E-20, 2.7105054312137611E-20,
    3.3881317890172014E-20, 3.3881317890172014E-20, 3.0493186101154812E-20,
    3.7269449679189215E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 3.7269449679189215E-20,
    4.4045713257223618E-20, 4.7433845046240819E-20, 4.7433845046240819E-20,
    5.082197683525802E-20, 5.7598240413292423E-20, 5.7598240413292423E-20,
    6.4374503991326826E-20, 6.4374503991326826E-20, 6.4374503991326826E-20,
    6.7762635780344027E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
    -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.10000000000000005,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000007,
    -3.0493186101154812E-20, 0.1000000000000001, -2.7105054312137611E-20,
    0.1000000000000001, -3.3881317890172014E-20, 0.10000000000000012,
    -3.3881317890172014E-20, 0.10000000000000012, -3.0493186101154812E-20,
    0.10000000000000013, -3.7269449679189215E-20, 0.10000000000000013,
    -4.0657581468206416E-20, 0.10000000000000014, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000017, -3.7269449679189215E-20,
    0.10000000000000019, -4.4045713257223618E-20, 0.10000000000000019,
    -4.7433845046240819E-20, 0.10000000000000021, -4.7433845046240819E-20,
    0.10000000000000021, -5.082197683525802E-20, 0.10000000000000023,
    -5.7598240413292423E-20, 0.10000000000000024, -5.7598240413292423E-20,
    0.10000000000000026, -6.4374503991326826E-20, 0.10000000000000026,
    -6.4374503991326826E-20, 0.10000000000000026, -6.4374503991326826E-20,
    0.10000000000000028, -6.7762635780344027E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 3.3881317890172014E-20,
    3.3881317890172014E-20, 3.0493186101154812E-20, 3.7269449679189215E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 3.7269449679189215E-20, 4.4045713257223618E-20,
    4.7433845046240819E-20, 4.7433845046240819E-20, 5.082197683525802E-20,
    5.7598240413292423E-20, 5.7598240413292423E-20, 6.4374503991326826E-20,
    6.4374503991326826E-20, 6.4374503991326826E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.1, -1.6940658945086007E-20, 0.10000000000000002,
    -2.0328790734103208E-20, 0.10000000000000003, -2.0328790734103208E-20,
    0.10000000000000005, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000007, -3.0493186101154812E-20, 0.1000000000000001,
    -2.7105054312137611E-20, 0.1000000000000001, -3.3881317890172014E-20,
    0.10000000000000012, -3.3881317890172014E-20, 0.10000000000000012,
    -3.0493186101154812E-20, 0.10000000000000013, -3.7269449679189215E-20,
    0.10000000000000013, -4.0657581468206416E-20, 0.10000000000000014,
    -4.0657581468206416E-20, 0.10000000000000016, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000017,
    -3.7269449679189215E-20, 0.10000000000000019, -4.4045713257223618E-20,
    0.10000000000000019, -4.7433845046240819E-20, 0.10000000000000021,
    -4.7433845046240819E-20, 0.10000000000000021, -5.082197683525802E-20,
    0.10000000000000023, -5.7598240413292423E-20, 0.10000000000000024,
    -5.7598240413292423E-20, 0.10000000000000026, -6.4374503991326826E-20,
    0.10000000000000026, -6.4374503991326826E-20, 0.10000000000000026,
    -6.4374503991326826E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 3.3881317890172014E-20,
    3.3881317890172014E-20, 3.0493186101154812E-20, 3.7269449679189215E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 3.7269449679189215E-20, 4.4045713257223618E-20,
    4.7433845046240819E-20, 4.7433845046240819E-20, 5.082197683525802E-20,
    5.7598240413292423E-20, 5.7598240413292423E-20, 6.4374503991326826E-20,
    6.4374503991326826E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.1, -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.10000000000000005,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000007,
    -3.0493186101154812E-20, 0.1000000000000001, -2.7105054312137611E-20,
    0.1000000000000001, -3.3881317890172014E-20, 0.10000000000000012,
    -3.3881317890172014E-20, 0.10000000000000012, -3.0493186101154812E-20,
    0.10000000000000013, -3.7269449679189215E-20, 0.10000000000000013,
    -4.0657581468206416E-20, 0.10000000000000014, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000017, -3.7269449679189215E-20,
    0.10000000000000019, -4.4045713257223618E-20, 0.10000000000000019,
    -4.7433845046240819E-20, 0.10000000000000021, -4.7433845046240819E-20,
    0.10000000000000021, -5.082197683525802E-20, 0.10000000000000023,
    -5.7598240413292423E-20, 0.10000000000000024, -5.7598240413292423E-20,
    0.10000000000000026, -6.4374503991326826E-20, 0.10000000000000026,
    -6.4374503991326826E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 3.3881317890172014E-20,
    3.3881317890172014E-20, 3.0493186101154812E-20, 3.7269449679189215E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 3.7269449679189215E-20, 4.4045713257223618E-20,
    4.7433845046240819E-20, 4.7433845046240819E-20, 5.082197683525802E-20,
    5.7598240413292423E-20, 5.7598240413292423E-20, 6.4374503991326826E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
    -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.10000000000000005,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000007,
    -3.0493186101154812E-20, 0.1000000000000001, -2.7105054312137611E-20,
    0.1000000000000001, -3.3881317890172014E-20, 0.10000000000000012,
    -3.3881317890172014E-20, 0.10000000000000012, -3.0493186101154812E-20,
    0.10000000000000013, -3.7269449679189215E-20, 0.10000000000000013,
    -4.0657581468206416E-20, 0.10000000000000014, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000017, -3.7269449679189215E-20,
    0.10000000000000019, -4.4045713257223618E-20, 0.10000000000000019,
    -4.7433845046240819E-20, 0.10000000000000021, -4.7433845046240819E-20,
    0.10000000000000021, -5.082197683525802E-20, 0.10000000000000023,
    -5.7598240413292423E-20, 0.10000000000000024, -5.7598240413292423E-20,
    0.10000000000000026, -6.4374503991326826E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20,
    2.0328790734103208E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 3.0493186101154812E-20, 2.7105054312137611E-20,
    3.3881317890172014E-20, 3.3881317890172014E-20, 3.0493186101154812E-20,
    3.7269449679189215E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 3.7269449679189215E-20,
    4.4045713257223618E-20, 4.7433845046240819E-20, 4.7433845046240819E-20,
    5.082197683525802E-20, 5.7598240413292423E-20, 5.7598240413292423E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
    -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.10000000000000005,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000007,
    -3.0493186101154812E-20, 0.1000000000000001, -2.7105054312137611E-20,
    0.1000000000000001, -3.3881317890172014E-20, 0.10000000000000012,
    -3.3881317890172014E-20, 0.10000000000000012, -3.0493186101154812E-20,
    0.10000000000000013, -3.7269449679189215E-20, 0.10000000000000013,
    -4.0657581468206416E-20, 0.10000000000000014, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000017, -3.7269449679189215E-20,
    0.10000000000000019, -4.4045713257223618E-20, 0.10000000000000019,
    -4.7433845046240819E-20, 0.10000000000000021, -4.7433845046240819E-20,
    0.10000000000000021, -5.082197683525802E-20, 0.10000000000000023,
    -5.7598240413292423E-20, 0.10000000000000024, -5.7598240413292423E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 3.3881317890172014E-20,
    3.3881317890172014E-20, 3.0493186101154812E-20, 3.7269449679189215E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 3.7269449679189215E-20, 4.4045713257223618E-20,
    4.7433845046240819E-20, 4.7433845046240819E-20, 5.082197683525802E-20,
    5.7598240413292423E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000007, -3.0493186101154812E-20,
    0.1000000000000001, -2.7105054312137611E-20, 0.1000000000000001,
    -3.3881317890172014E-20, 0.10000000000000012, -3.3881317890172014E-20,
    0.10000000000000012, -3.0493186101154812E-20, 0.10000000000000013,
    -3.7269449679189215E-20, 0.10000000000000013, -4.0657581468206416E-20,
    0.10000000000000014, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000016, -4.0657581468206416E-20,
    0.10000000000000017, -3.7269449679189215E-20, 0.10000000000000019,
    -4.4045713257223618E-20, 0.10000000000000019, -4.7433845046240819E-20,
    0.10000000000000021, -4.7433845046240819E-20, 0.10000000000000021,
    -5.082197683525802E-20, 0.10000000000000023, -5.7598240413292423E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 3.3881317890172014E-20,
    3.3881317890172014E-20, 3.0493186101154812E-20, 3.7269449679189215E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 3.7269449679189215E-20, 4.4045713257223618E-20,
    4.7433845046240819E-20, 4.7433845046240819E-20, 5.082197683525802E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.1, -1.6940658945086007E-20, 0.10000000000000002,
    -2.0328790734103208E-20, 0.10000000000000003, -2.0328790734103208E-20,
    0.10000000000000005, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000007, -3.0493186101154812E-20, 0.1000000000000001,
    -2.7105054312137611E-20, 0.1000000000000001, -3.3881317890172014E-20,
    0.10000000000000012, -3.3881317890172014E-20, 0.10000000000000012,
    -3.0493186101154812E-20, 0.10000000000000013, -3.7269449679189215E-20,
    0.10000000000000013, -4.0657581468206416E-20, 0.10000000000000014,
    -4.0657581468206416E-20, 0.10000000000000016, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000017,
    -3.7269449679189215E-20, 0.10000000000000019, -4.4045713257223618E-20,
    0.10000000000000019, -4.7433845046240819E-20, 0.10000000000000021,
    -4.7433845046240819E-20, 0.10000000000000021, -5.082197683525802E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 3.3881317890172014E-20,
    3.3881317890172014E-20, 3.0493186101154812E-20, 3.7269449679189215E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 3.7269449679189215E-20, 4.4045713257223618E-20,
    4.7433845046240819E-20, 4.7433845046240819E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
    -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.10000000000000005,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000007,
    -3.0493186101154812E-20, 0.1000000000000001, -2.7105054312137611E-20,
    0.1000000000000001, -3.3881317890172014E-20, 0.10000000000000012,
    -3.3881317890172014E-20, 0.10000000000000012, -3.0493186101154812E-20,
    0.10000000000000013, -3.7269449679189215E-20, 0.10000000000000013,
    -4.0657581468206416E-20, 0.10000000000000014, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000017, -3.7269449679189215E-20,
    0.10000000000000019, -4.4045713257223618E-20, 0.10000000000000019,
    -4.7433845046240819E-20, 0.10000000000000021, -4.7433845046240819E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 3.3881317890172014E-20,
    3.3881317890172014E-20, 3.0493186101154812E-20, 3.7269449679189215E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 3.7269449679189215E-20, 4.4045713257223618E-20,
    4.7433845046240819E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
    -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.10000000000000005,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000007,
    -3.0493186101154812E-20, 0.1000000000000001, -2.7105054312137611E-20,
    0.1000000000000001, -3.3881317890172014E-20, 0.10000000000000012,
    -3.3881317890172014E-20, 0.10000000000000012, -3.0493186101154812E-20,
    0.10000000000000013, -3.7269449679189215E-20, 0.10000000000000013,
    -4.0657581468206416E-20, 0.10000000000000014, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000017, -3.7269449679189215E-20,
    0.10000000000000019, -4.4045713257223618E-20, 0.10000000000000019,
    -4.7433845046240819E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20,
    2.0328790734103208E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 3.0493186101154812E-20, 2.7105054312137611E-20,
    3.3881317890172014E-20, 3.3881317890172014E-20, 3.0493186101154812E-20,
    3.7269449679189215E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 3.7269449679189215E-20,
    4.4045713257223618E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
    -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.10000000000000005,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000007,
    -3.0493186101154812E-20, 0.1000000000000001, -2.7105054312137611E-20,
    0.1000000000000001, -3.3881317890172014E-20, 0.10000000000000012,
    -3.3881317890172014E-20, 0.10000000000000012, -3.0493186101154812E-20,
    0.10000000000000013, -3.7269449679189215E-20, 0.10000000000000013,
    -4.0657581468206416E-20, 0.10000000000000014, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.10000000000000017, -3.7269449679189215E-20,
    0.10000000000000019, -4.4045713257223618E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 3.3881317890172014E-20,
    3.3881317890172014E-20, 3.0493186101154812E-20, 3.7269449679189215E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 3.7269449679189215E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20, 0.10000000000000002,
    -2.0328790734103208E-20, 0.10000000000000003, -2.0328790734103208E-20,
    0.10000000000000005, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000007, -3.0493186101154812E-20, 0.1000000000000001,
    -2.7105054312137611E-20, 0.1000000000000001, -3.3881317890172014E-20,
    0.10000000000000012, -3.3881317890172014E-20, 0.10000000000000012,
    -3.0493186101154812E-20, 0.10000000000000013, -3.7269449679189215E-20,
    0.10000000000000013, -4.0657581468206416E-20, 0.10000000000000014,
    -4.0657581468206416E-20, 0.10000000000000016, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000017,
    -3.7269449679189215E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20,
    2.0328790734103208E-20, 2.0328790734103208E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 3.0493186101154812E-20,
    2.7105054312137611E-20, 3.3881317890172014E-20, 3.3881317890172014E-20,
    3.0493186101154812E-20, 3.7269449679189215E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 4.0657581468206416E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
    -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.10000000000000005,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000007,
    -3.0493186101154812E-20, 0.1000000000000001, -2.7105054312137611E-20,
    0.1000000000000001, -3.3881317890172014E-20, 0.10000000000000012,
    -3.3881317890172014E-20, 0.10000000000000012, -3.0493186101154812E-20,
    0.10000000000000013, -3.7269449679189215E-20, 0.10000000000000013,
    -4.0657581468206416E-20, 0.10000000000000014, -4.0657581468206416E-20,
    0.10000000000000016, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20,
    2.0328790734103208E-20, 2.0328790734103208E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 3.0493186101154812E-20,
    2.7105054312137611E-20, 3.3881317890172014E-20, 3.3881317890172014E-20,
    3.0493186101154812E-20, 3.7269449679189215E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000007, -3.0493186101154812E-20,
    0.1000000000000001, -2.7105054312137611E-20, 0.1000000000000001,
    -3.3881317890172014E-20, 0.10000000000000012, -3.3881317890172014E-20,
    0.10000000000000012, -3.0493186101154812E-20, 0.10000000000000013,
    -3.7269449679189215E-20, 0.10000000000000013, -4.0657581468206416E-20,
    0.10000000000000014, -4.0657581468206416E-20, 0.10000000000000016,
    -4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20,
    2.0328790734103208E-20, 2.0328790734103208E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 3.0493186101154812E-20,
    2.7105054312137611E-20, 3.3881317890172014E-20, 3.3881317890172014E-20,
    3.0493186101154812E-20, 3.7269449679189215E-20, 4.0657581468206416E-20,
    4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000007, -3.0493186101154812E-20,
    0.1000000000000001, -2.7105054312137611E-20, 0.1000000000000001,
    -3.3881317890172014E-20, 0.10000000000000012, -3.3881317890172014E-20,
    0.10000000000000012, -3.0493186101154812E-20, 0.10000000000000013,
    -3.7269449679189215E-20, 0.10000000000000013, -4.0657581468206416E-20,
    0.10000000000000014, -4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 3.3881317890172014E-20,
    3.3881317890172014E-20, 3.0493186101154812E-20, 3.7269449679189215E-20,
    4.0657581468206416E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000007, -3.0493186101154812E-20,
    0.1000000000000001, -2.7105054312137611E-20, 0.1000000000000001,
    -3.3881317890172014E-20, 0.10000000000000012, -3.3881317890172014E-20,
    0.10000000000000012, -3.0493186101154812E-20, 0.10000000000000013,
    -3.7269449679189215E-20, 0.10000000000000013, -4.0657581468206416E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20,
    2.0328790734103208E-20, 2.0328790734103208E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 3.0493186101154812E-20,
    2.7105054312137611E-20, 3.3881317890172014E-20, 3.3881317890172014E-20,
    3.0493186101154812E-20, 3.7269449679189215E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.1, -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.10000000000000005,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000007,
    -3.0493186101154812E-20, 0.1000000000000001, -2.7105054312137611E-20,
    0.1000000000000001, -3.3881317890172014E-20, 0.10000000000000012,
    -3.3881317890172014E-20, 0.10000000000000012, -3.0493186101154812E-20,
    0.10000000000000013, -3.7269449679189215E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20,
    2.0328790734103208E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 3.0493186101154812E-20, 2.7105054312137611E-20,
    3.3881317890172014E-20, 3.3881317890172014E-20, 3.0493186101154812E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000007, -3.0493186101154812E-20,
    0.1000000000000001, -2.7105054312137611E-20, 0.1000000000000001,
    -3.3881317890172014E-20, 0.10000000000000012, -3.3881317890172014E-20,
    0.10000000000000012, -3.0493186101154812E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20,
    2.0328790734103208E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 3.0493186101154812E-20, 2.7105054312137611E-20,
    3.3881317890172014E-20, 3.3881317890172014E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20, 0.10000000000000002,
    -2.0328790734103208E-20, 0.10000000000000003, -2.0328790734103208E-20,
    0.10000000000000005, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000007, -3.0493186101154812E-20, 0.1000000000000001,
    -2.7105054312137611E-20, 0.1000000000000001, -3.3881317890172014E-20,
    0.10000000000000012, -3.3881317890172014E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20,
    2.0328790734103208E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 3.0493186101154812E-20, 2.7105054312137611E-20,
    3.3881317890172014E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.1, -1.6940658945086007E-20, 0.10000000000000002,
    -2.0328790734103208E-20, 0.10000000000000003, -2.0328790734103208E-20,
    0.10000000000000005, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000007, -3.0493186101154812E-20, 0.1000000000000001,
    -2.7105054312137611E-20, 0.1000000000000001, -3.3881317890172014E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    1.6940658945086007E-20, 2.0328790734103208E-20, 2.0328790734103208E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    3.0493186101154812E-20, 2.7105054312137611E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.10000000000000007, -3.0493186101154812E-20,
    0.1000000000000001, -2.7105054312137611E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20,
    2.0328790734103208E-20, 2.0328790734103208E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 2.3716922523120409E-20, 3.0493186101154812E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.1, -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.10000000000000005,
    -2.3716922523120409E-20, 0.10000000000000006, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000007,
    -3.0493186101154812E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20,
    2.0328790734103208E-20, 2.3716922523120409E-20, 2.3716922523120409E-20,
    2.3716922523120409E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20,
    0.10000000000000006, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20,
    2.0328790734103208E-20, 2.3716922523120409E-20, 2.3716922523120409E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20, 0.10000000000000002,
    -2.0328790734103208E-20, 0.10000000000000003, -2.0328790734103208E-20,
    0.10000000000000005, -2.3716922523120409E-20, 0.10000000000000006,
    -2.3716922523120409E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20,
    2.0328790734103208E-20, 2.0328790734103208E-20, 2.3716922523120409E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20,
    0.10000000000000002, -2.0328790734103208E-20, 0.10000000000000003,
    -2.0328790734103208E-20, 0.10000000000000005, -2.3716922523120409E-20, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, 1.6940658945086007E-20, 2.0328790734103208E-20,
    2.0328790734103208E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
    -1.6940658945086007E-20, 0.10000000000000002, -2.0328790734103208E-20,
    0.10000000000000003, -2.0328790734103208E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    1.6940658945086007E-20, 2.0328790734103208E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.1, -1.6940658945086007E-20, 0.10000000000000002,
    -2.0328790734103208E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.6940658945086007E-20,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1,
    -1.6940658945086007E-20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Mx[384] = { -0.0051781079403026339,
    -0.0042394762134830462, -0.0034709895529211579, -0.0028418058905889396,
    -0.00232667387690332, -0.0019049194554039152, -0.0015596161402757387,
    -0.0012769056970405234, -0.0010454419629475453, -0.00085593548562335961,
    -0.00070078070473057351, -0.00057375071412657346, -0.00046974735425587161,
    -0.00038459660510629148, -0.00031488106812990144, -0.00025780281403998566,
    -0.00021107109208457254, -0.00017281039417538628, -0.00014148518416290908,
    -0.00011583827137906764, -9.4840355161427835E-5, -7.7648715403491234E-5,
    -6.35733912378304E-5, -5.20494904838622E-5, -4.2614518541169891E-5,
    -3.4889816857259789E-5, -2.8565366030289315E-5, -2.3387343641919111E-5,
    -1.9147937472434125E-5, -1.5677005266688235E-5, -0.00057322369001704034,
    0.0051781079403026339, -0.0010425395534268307, 0.0042394762134830462,
    -0.0014267828837077712, 0.0034709895529211579, -0.0017413747148738762,
    0.0028418058905889396, -0.0019989407217166823, 0.00232667387690332,
    -0.0022098179324663812, 0.0019049194554039152, -0.0023824695900304654,
    0.0015596161402757387, -0.0025238248116480691, 0.0012769056970405234,
    -0.002639556678694554, 0.0010454419629475453, -0.0027343099173566428,
    0.00085593548562335961, -0.0028118873078030321, 0.00070078070473057351,
    -0.002875402303105028, 0.00057375071412657346, -0.0029274039830403751,
    0.00046974735425587161, -0.0029699793576151612, 0.00038459660510629148,
    -0.0030048371261033519, 0.00031488106812990144, -0.0030333762531483061,
    0.00025780281403998566, -0.0030567421141260084, 0.00021107109208457254,
    -0.0030758724630805976, 0.00017281039417538628, -0.0030915350680868325,
    0.00014148518416290908, -0.0031043585244787493, 0.00011583827137906764,
    -0.0031148574825875649, 9.4840355161427835E-5, -0.0031234533024665292,
    7.7648715403491234E-5, -0.0031304909645493558, 6.35733912378304E-5,
    -0.0031362529149263356, 5.20494904838622E-5, -0.0031409704008976778,
    4.2614518541169891E-5, -0.0031448327517396282, 3.4889816857259789E-5,
    -0.0031479949771531096, 2.8565366030289315E-5, -0.0031505839883472909,
    2.3387343641919111E-5, -0.0031527036914320294, 1.9147937472434125E-5,
    -0.0031544391575348981, 1.5677005266688235E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.013222333301565128, 0.013691649128446955, 0.014075892428821329,
    0.014390484235502009, 0.014648050222297844, 0.014858927416634473,
    0.015031579060760672, 0.015172934271376268, 0.01528866612941507,
    0.015383419360702295, 0.015460996745110655, 0.015524511735469131,
    0.015576513411357068, 0.015619088782618116, 0.015653946548393249,
    0.01568248567321694, 0.015705851532376025, 0.015724981879841657,
    0.015740644483628836, 0.015753467939022676, 0.015763966896314336,
    0.015772562715524272, 0.015779600377059343, 0.015785362326987861,
    0.015790079812592034, 0.015793942163133373, 0.015797104388300735,
    0.015799693399293411, 0.01580181310221317, 0.015803548568180968,
    0.0012833675034005457, -0.013222333301565128, 0.0026298482969800297,
    -0.013691649128446955, 0.0040188653538132408, -0.014075892428821329,
    0.005442708157493298, -0.014390484235502009, 0.0068950638711157784,
    -0.014648050222297844, 0.0083707639809678623, -0.014858927416634473,
    0.0098655768659251623, -0.015031579060760672, 0.011376037967637764,
    -0.015172934271376268, 0.012899310745638761, -0.01528866612941507,
    0.014433072837015549, -0.015383419360702295, 0.015975422851831768,
    -0.015460996745110655, 0.01752480406367293, -0.015524511735469131,
    0.01907994193274936, -0.015576513411357068, 0.020639792954139233,
    -0.015619088782618116, 0.022203502778272068, -0.015653946548393249,
    0.023770371922880625, -0.01568248567321694, 0.02533982770032148,
    -0.015705851532376025, 0.026911401233609064, -0.015724981879841657,
    0.028484708638735871, -0.015740644483628836, 0.030059435618059353,
    -0.015753467939022676, 0.031635324846433922, -0.015763966896314336,
    0.033212165643849351, -0.015772562715524272, 0.034789785520102215,
    -0.015779600377059343, 0.036368043252158359, -0.015785362326987861,
    0.037946823216376671, -0.015790079812592034, 0.039526030748126405,
    -0.015793942163133373, 0.041105588342563122, -0.015797104388300735,
    0.04268543254408718, -0.015799693399293411, 0.044265511399647704,
    -0.01580181310221317, 0.0458457823736842, -0.015803548568180968, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 5.2167886576861573E-6, 5.4019542730666079E-6,
    5.5535550567911924E-6, 5.6776752806172337E-6, 5.7792963249425364E-6,
    5.8624965990915628E-6, 5.9306152222018894E-6, 5.9863860337996428E-6,
    6.032047312378842E-6, 6.0694316053764926E-6, 6.1000392757357475E-6,
    6.1250987167389431E-6, 6.1456156517432022E-6, 6.1624134973900915E-6,
    6.1761664102066576E-6, 6.18742634287398E-6, 6.1966451960263048E-6,
    6.2041929546102232E-6, 6.2103725366796848E-6, 6.2154319505611228E-6,
    6.2195742482984058E-6, 6.2229656748443252E-6, 6.2257423400542747E-6,
    6.2280156812526614E-6, 6.22987693560402E-6, 6.231400801780778E-6,
    6.2326484378832653E-6, 6.2336699159290219E-6, 6.234506231418677E-6,
    6.2351909486293334E-6, 0.028284779793491491, -5.2167886576861573E-6,
    0.028285311039041847, -5.4019542730666079E-6, 0.028285859067008031,
    -5.5535550567911924E-6, 0.028286420835254164, -5.6776752806172337E-6,
    0.028286993853090046, -5.7792963249425364E-6, 0.028287576081311019,
    -5.8624965990915628E-6, 0.028288165850357511, -5.9306152222018894E-6,
    0.028288761793309757, -5.9863860337996428E-6, 0.028289362791028512,
    -6.032047312378842E-6, 0.028289967927240058, -6.0694316053764926E-6,
    0.028290576451762924, -6.1000392757357475E-6, 0.028291187750400466,
    -6.1250987167389431E-6, 0.028291801320291011, -6.1456156517432022E-6,
    0.028292416749726233, -6.1624134973900915E-6, 0.028293033701627869,
    -6.1761664102066576E-6, 0.028293651900019583, -6.18742634287398E-6,
    0.028294271118951053, -6.1966451960263048E-6, 0.028294891173429804,
    -6.2041929546102232E-6, 0.028295511911996817, -6.2103725366796848E-6,
    0.028296133210647924, -6.2154319505611228E-6, 0.0282967549678571,
    -6.2195742482984058E-6, 0.028297377100501876, -6.2229656748443252E-6,
    0.028297999540527315, -6.2257423400542747E-6, 0.02829862223221476,
    -6.2280156812526614E-6, 0.028299245129945631, -6.22987693560402E-6,
    0.028299868196370589, -6.231400801780778E-6, 0.028300491400910583,
    -6.2326484378832653E-6, 0.028301114718529606, -6.2336699159290219E-6,
    0.028301738128729905, -6.234506231418677E-6, 0.028302361614729345,
    -6.2351909486293334E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Mu1[96] = { -0.00093653765389917893,
    -0.0035160023017821343, -0.0074405818047015727, -0.012466448205861413,
    -0.01839397205857253, -0.025059710595610598, -0.032329848197080895,
    -0.04009482589973342, -0.048264944411080056, -0.056766764161831441,
    -0.065540157918117584, -0.074535897664471593, -0.083713678910717737,
    -0.093040503131262026, -0.10248935341839441, -0.1120381101989196,
    -0.12166866349801768, -0.13136618612236609, -0.14111853859280982,
    -0.15091578194443833, -0.16074977884102559, -0.17061386699515521,
    -0.18050259178723355, -0.19041148735245295, -0.2003368973499563,
    -0.21027582822104013, -0.22022582904713281, -0.23018489318582641,
    -0.24015137773727113, -0.25012393760883572, -3.1731173050448962E-5,
    0.00093653765389917893, -0.000241998849108999, 0.0035160023017821343,
    -0.00077970909764929713, 0.0074405818047015727, -0.0017667758970693849,
    0.012466448205861413, -0.003303013970713824, 0.01839397205857253,
    -0.0054701447021947777, 0.025059710595610598, -0.0083350759014596074,
    0.032329848197080895, -0.011952587050133315, 0.04009482589973342,
    -0.016367527794459956, 0.048264944411080056, -0.021616617919084216,
    0.056766764161831441, -0.027729921040941088, 0.065540157918117584,
    -0.034732051167764018, 0.074535897664471593, -0.042643160544640867,
    0.083713678910717737, -0.051479748434368644, 0.093040503131262026,
    -0.061255323290802362, 0.10248935341839441, -0.071980944900539662,
    0.1120381101989196, -0.083665668250990519, 0.12166866349801768,
    -0.0963169069388162, 0.13136618612236609, -0.10994073070359421,
    0.14111853859280982, -0.12454210902777982, 0.15091578194443833,
    -0.14012511057948604, 0.16074977884102559, -0.1566930665024211,
    0.17061386699515521, -0.17424870410638177, 0.18050259178723355,
    -0.1927942563237719, 0.19041148735245295, -0.21233155132502005,
    0.2003368973499563, -0.23286208588947793, 0.21027582822104013,
    -0.25438708547643141, 0.22022582904713281, -0.27690755340708439,
    0.23018489318582641, -0.30042431113136181, 0.24015137773727113,
    -0.32493803119557929, 0.25012393760883572, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0 };

  static const uint8_T b_Mrows_0[96] = { 2U, 4U, 6U, 8U, 10U, 12U, 14U, 16U, 18U,
    20U, 22U, 24U, 26U, 28U, 30U, 32U, 34U, 36U, 38U, 40U, 42U, 44U, 46U, 48U,
    50U, 52U, 54U, 56U, 58U, 60U, 61U, 62U, 63U, 64U, 65U, 66U, 67U, 68U, 69U,
    70U, 71U, 72U, 73U, 74U, 75U, 76U, 77U, 78U, 79U, 80U, 81U, 82U, 83U, 84U,
    85U, 86U, 87U, 88U, 89U, 90U, 91U, 92U, 93U, 94U, 95U, 96U, 97U, 98U, 99U,
    100U, 101U, 102U, 103U, 104U, 105U, 106U, 107U, 108U, 109U, 110U, 111U, 112U,
    113U, 114U, 115U, 116U, 117U, 118U, 119U, 120U, 121U, 122U, 123U, 151U, 152U,
    153U };

  static const real_T b_Linv[16] = { 7.9998080541108605, -2.9530049756913419,
    -2.1896728723561925, 0.0, 0.0, 8.7034491930031113, -2.0693880092142352, 0.0,
    0.0, 0.0, 9.0725349962825419, 0.0, 0.0, 0.0, 0.0, 0.001 };

  static const real_T b_Hinv[16] = { 77.51183457700742, -21.170045986459385,
    -19.865883764862073, 0.0, -21.170045986459385, 80.03239458786615,
    -18.774595134483608, 0.0, -19.865883764862073, -18.774595134483608,
    82.310891258771463, 0.0, 0.0, 0.0, 0.0, 1.0E-6 };

  static const real_T b_Ac[384] = { -0.00093653765389917893,
    -0.0035160023017821343, -0.0074405818047015727, -0.012466448205861413,
    -0.01839397205857253, -0.025059710595610598, -0.032329848197080895,
    -0.04009482589973342, -0.048264944411080056, -0.056766764161831441,
    -0.065540157918117584, -0.074535897664471593, -0.083713678910717737,
    -0.093040503131262026, -0.10248935341839441, -0.1120381101989196,
    -0.12166866349801768, -0.13136618612236609, -0.14111853859280982,
    -0.15091578194443833, -0.16074977884102559, -0.17061386699515521,
    -0.18050259178723355, -0.19041148735245295, -0.2003368973499563,
    -0.21027582822104013, -0.22022582904713281, -0.23018489318582641,
    -0.24015137773727113, -0.25012393760883572, -3.1731173050448962E-5,
    0.00093653765389917893, -0.000241998849108999, 0.0035160023017821343,
    -0.00077970909764929713, 0.0074405818047015727, -0.0017667758970693849,
    0.012466448205861413, -0.003303013970713824, 0.01839397205857253,
    -0.0054701447021947777, 0.025059710595610598, -0.0083350759014596074,
    0.032329848197080895, -0.011952587050133315, 0.04009482589973342,
    -0.016367527794459956, 0.048264944411080056, -0.021616617919084216,
    0.056766764161831441, -0.027729921040941088, 0.065540157918117584,
    -0.034732051167764018, 0.074535897664471593, -0.042643160544640867,
    0.083713678910717737, -0.051479748434368644, 0.093040503131262026,
    -0.061255323290802362, 0.10248935341839441, -0.071980944900539662,
    0.1120381101989196, -0.083665668250990519, 0.12166866349801768,
    -0.0963169069388162, 0.13136618612236609, -0.10994073070359421,
    0.14111853859280982, -0.12454210902777982, 0.15091578194443833,
    -0.14012511057948604, 0.16074977884102559, -0.1566930665024211,
    0.17061386699515521, -0.17424870410638177, 0.18050259178723355,
    -0.1927942563237719, 0.19041148735245295, -0.21233155132502005,
    0.2003368973499563, -0.23286208588947793, 0.21027582822104013,
    -0.25438708547643141, 0.22022582904713281, -0.27690755340708439,
    0.23018489318582641, -0.30042431113136181, 0.24015137773727113,
    -0.32493803119557929, 0.25012393760883572, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0,
    -0.0, -0.00093653765389917893, -0.0035160023017821343,
    -0.0074405818047015727, -0.012466448205861413, -0.01839397205857253,
    -0.025059710595610598, -0.032329848197080895, -0.04009482589973342,
    -0.048264944411080056, -0.056766764161831441, -0.065540157918117584,
    -0.074535897664471593, -0.083713678910717737, -0.093040503131262026,
    -0.10248935341839441, -0.1120381101989196, -0.12166866349801768,
    -0.13136618612236609, -0.14111853859280982, -0.15091578194443833,
    -0.16074977884102559, -0.17061386699515521, -0.18050259178723355,
    -0.19041148735245295, -0.2003368973499563, -0.21027582822104013,
    -0.22022582904713281, -0.23018489318582641, -0.24015137773727113, 0.0, 0.0,
    -3.1731173050448962E-5, 0.00093653765389917893, -0.000241998849108999,
    0.0035160023017821343, -0.00077970909764929713, 0.0074405818047015727,
    -0.0017667758970693849, 0.012466448205861413, -0.003303013970713824,
    0.01839397205857253, -0.0054701447021947777, 0.025059710595610598,
    -0.0083350759014596074, 0.032329848197080895, -0.011952587050133315,
    0.04009482589973342, -0.016367527794459956, 0.048264944411080056,
    -0.021616617919084216, 0.056766764161831441, -0.027729921040941088,
    0.065540157918117584, -0.034732051167764018, 0.074535897664471593,
    -0.042643160544640867, 0.083713678910717737, -0.051479748434368644,
    0.093040503131262026, -0.061255323290802362, 0.10248935341839441,
    -0.071980944900539662, 0.1120381101989196, -0.083665668250990519,
    0.12166866349801768, -0.0963169069388162, 0.13136618612236609,
    -0.10994073070359421, 0.14111853859280982, -0.12454210902777982,
    0.15091578194443833, -0.14012511057948604, 0.16074977884102559,
    -0.1566930665024211, 0.17061386699515521, -0.17424870410638177,
    0.18050259178723355, -0.1927942563237719, 0.19041148735245295,
    -0.21233155132502005, 0.2003368973499563, -0.23286208588947793,
    0.21027582822104013, -0.25438708547643141, 0.22022582904713281,
    -0.27690755340708439, 0.23018489318582641, -0.30042431113136181,
    0.24015137773727113, -0.0, -1.0, -1.0, 0.0, 1.0, 1.0, -0.0, -0.0,
    -0.00093653765389917893, -0.0035160023017821343, -0.0074405818047015727,
    -0.012466448205861413, -0.01839397205857253, -0.025059710595610598,
    -0.032329848197080895, -0.04009482589973342, -0.048264944411080056,
    -0.056766764161831441, -0.065540157918117584, -0.074535897664471593,
    -0.083713678910717737, -0.093040503131262026, -0.10248935341839441,
    -0.1120381101989196, -0.12166866349801768, -0.13136618612236609,
    -0.14111853859280982, -0.15091578194443833, -0.16074977884102559,
    -0.17061386699515521, -0.18050259178723355, -0.19041148735245295,
    -0.2003368973499563, -0.21027582822104013, -0.22022582904713281,
    -0.23018489318582641, 0.0, 0.0, 0.0, 0.0, -3.1731173050448962E-5,
    0.00093653765389917893, -0.000241998849108999, 0.0035160023017821343,
    -0.00077970909764929713, 0.0074405818047015727, -0.0017667758970693849,
    0.012466448205861413, -0.003303013970713824, 0.01839397205857253,
    -0.0054701447021947777, 0.025059710595610598, -0.0083350759014596074,
    0.032329848197080895, -0.011952587050133315, 0.04009482589973342,
    -0.016367527794459956, 0.048264944411080056, -0.021616617919084216,
    0.056766764161831441, -0.027729921040941088, 0.065540157918117584,
    -0.034732051167764018, 0.074535897664471593, -0.042643160544640867,
    0.083713678910717737, -0.051479748434368644, 0.093040503131262026,
    -0.061255323290802362, 0.10248935341839441, -0.071980944900539662,
    0.1120381101989196, -0.083665668250990519, 0.12166866349801768,
    -0.0963169069388162, 0.13136618612236609, -0.10994073070359421,
    0.14111853859280982, -0.12454210902777982, 0.15091578194443833,
    -0.14012511057948604, 0.16074977884102559, -0.1566930665024211,
    0.17061386699515521, -0.17424870410638177, 0.18050259178723355,
    -0.1927942563237719, 0.19041148735245295, -0.21233155132502005,
    0.2003368973499563, -0.23286208588947793, 0.21027582822104013,
    -0.25438708547643141, 0.22022582904713281, -0.27690755340708439,
    0.23018489318582641, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0 };

  static const real_T b_Kr[180] = { 0.0, -9.36537653899179E-6, 0.0,
    -3.5160023017821348E-5, 0.0, -7.4405818047015737E-5, 0.0,
    -0.00012466448205861416, 0.0, -0.00018393972058572533, 0.0,
    -0.000250597105956106, 0.0, -0.000323298481970809, 0.0,
    -0.00040094825899733431, 0.0, -0.00048264944411080063, 0.0,
    -0.00056766764161831457, 0.0, -0.000655401579181176, 0.0,
    -0.00074535897664471606, 0.0, -0.00083713678910717754, 0.0,
    -0.00093040503131262039, 0.0, -0.0010248935341839443, 0.0,
    -0.0011203811019891961, 0.0, -0.001216686634980177, 0.0,
    -0.001313661861223661, 0.0, -0.0014111853859280984, 0.0,
    -0.0015091578194443836, 0.0, -0.0016074977884102562, 0.0,
    -0.0017061386699515524, 0.0, -0.0018050259178723359, 0.0,
    -0.0019041148735245298, 0.0, -0.0020033689734995632, 0.0,
    -0.0021027582822104019, 0.0, -0.0022022582904713283, 0.0,
    -0.0023018489318582646, 0.0, -0.0024015137773727119, 0.0,
    -0.0025012393760883577, -0.0, -0.0, 0.0, -9.36537653899179E-6, 0.0,
    -3.5160023017821348E-5, 0.0, -7.4405818047015737E-5, 0.0,
    -0.00012466448205861416, 0.0, -0.00018393972058572533, 0.0,
    -0.000250597105956106, 0.0, -0.000323298481970809, 0.0,
    -0.00040094825899733431, 0.0, -0.00048264944411080063, 0.0,
    -0.00056766764161831457, 0.0, -0.000655401579181176, 0.0,
    -0.00074535897664471606, 0.0, -0.00083713678910717754, 0.0,
    -0.00093040503131262039, 0.0, -0.0010248935341839443, 0.0,
    -0.0011203811019891961, 0.0, -0.001216686634980177, 0.0,
    -0.001313661861223661, 0.0, -0.0014111853859280984, 0.0,
    -0.0015091578194443836, 0.0, -0.0016074977884102562, 0.0,
    -0.0017061386699515524, 0.0, -0.0018050259178723359, 0.0,
    -0.0019041148735245298, 0.0, -0.0020033689734995632, 0.0,
    -0.0021027582822104019, 0.0, -0.0022022582904713283, 0.0,
    -0.0023018489318582646, 0.0, -0.0024015137773727119, -0.0, -0.0, -0.0, -0.0,
    0.0, -9.36537653899179E-6, 0.0, -3.5160023017821348E-5, 0.0,
    -7.4405818047015737E-5, 0.0, -0.00012466448205861416, 0.0,
    -0.00018393972058572533, 0.0, -0.000250597105956106, 0.0,
    -0.000323298481970809, 0.0, -0.00040094825899733431, 0.0,
    -0.00048264944411080063, 0.0, -0.00056766764161831457, 0.0,
    -0.000655401579181176, 0.0, -0.00074535897664471606, 0.0,
    -0.00083713678910717754, 0.0, -0.00093040503131262039, 0.0,
    -0.0010248935341839443, 0.0, -0.0011203811019891961, 0.0,
    -0.001216686634980177, 0.0, -0.001313661861223661, 0.0,
    -0.0014111853859280984, 0.0, -0.0015091578194443836, 0.0,
    -0.0016074977884102562, 0.0, -0.0017061386699515524, 0.0,
    -0.0018050259178723359, 0.0, -0.0019041148735245298, 0.0,
    -0.0020033689734995632, 0.0, -0.0021027582822104019, 0.0,
    -0.0022022582904713283, 0.0, -0.0023018489318582646 };

  static const real_T b_Kv[186] = { -1.853899189406451E-21, 0.0,
    -1.7798858199957204E-21, 0.0, -1.7091012039851541E-21, 0.0,
    -1.6415445464954563E-21, 0.0, -1.5684640325036183E-21, 0.0,
    -1.4984452197628245E-21, 0.0, -1.4228278469899485E-21, 0.0,
    -1.3414598491599051E-21, 0.0, -1.2718418872415314E-21, 0.0,
    -1.1969232476732372E-21, 0.0, -1.134245862615767E-21, 0.0,
    -1.0756053891422323E-21, 0.0, -1.0124718723758073E-21, 0.0,
    -9.5361251379423863E-22, 0.0, -9.08140930699146E-22, 0.0,
    -8.50265455527609E-22, 0.0, -7.8843589275488572E-22, 0.0,
    -7.2263199994526535E-22, 0.0, -6.5283720399411481E-22, 0.0,
    -5.8785051663166769E-22, 0.0, -5.3681159919151778E-22, 0.0,
    -4.7395010716065694E-22, 0.0, -4.0773330285105475E-22, 0.0,
    -3.557803983761415E-22, 0.0, -2.9232516100154664E-22, 0.0,
    -2.4346338076020297E-22, 0.0, -1.9224308336486977E-22, 0.0,
    -1.3866188057539606E-22, 0.0, -9.1530397706795409E-23, 0.0,
    -4.2372643210332577E-23, 0.0, 0.0, 0.0, -1.736303138489891E-21, 0.0,
    -1.6674595592809876E-21, 0.0, -1.6019207185123236E-21, 0.0,
    -1.5396106311438238E-21, 0.0, -1.472053973654126E-21, 0.0,
    -1.4074479883043546E-21, 0.0, -1.3374291755635608E-21, 0.0,
    -1.2618118027906848E-21, 0.0, -1.1973928622447743E-21, 0.0,
    -1.1277749003264006E-21, 0.0, -1.0698053180422395E-21, 0.0,
    -1.0156024616268357E-21, 0.0, -9.56961988153301E-22, 0.0,
    -9.0230300002894259E-22, 0.0, -8.60392698731507E-22, 0.0,
    -8.0644658699434782E-22, 0.0, -7.4857111182281084E-22, 0.0,
    -6.8674154905008755E-22, 0.0, -6.2093765624046718E-22, 0.0,
    -5.596173889313831E-22, 0.0, -5.11579758853069E-22, 0.0,
    -4.5206631277085264E-22, 0.0, -3.8920482073999175E-22, 0.0,
    -3.3993707371452263E-22, 0.0, -2.7950964059754287E-22, 0.0,
    -2.3300346050708104E-22, 0.0, -1.8414168026573734E-22, 0.0,
    -1.3292138287040417E-22, 0.0, -8.7814708722996963E-23, 0.0,
    -4.0683225854396318E-23, 0.0, 0.0, 0.0, -1.6212951513054073E-21, 0.0,
    -1.5572969447305471E-21, 0.0, -1.4965900106925231E-21, 0.0,
    -1.4391878150947382E-21, 0.0, -1.3768777277262384E-21, 0.0,
    -1.31745771540742E-21, 0.0, -1.2528517300576486E-21, 0.0,
    -1.1828329173168548E-21, 0.0, -1.1234888348857374E-21, 0.0,
    -1.0590698943398269E-21, 0.0, -1.0057252227632117E-21, 0.0,
    -9.5589228564992981E-22, 0.0, -9.01689429234526E-22, 0.0,
    -8.5118560093187054E-22, 0.0, -8.1279990314927064E-22, 0.0,
    -7.6275295668095589E-22, 0.0, -7.0880684494379669E-22, 0.0,
    -6.5093136977225971E-22, 0.0, -5.8910180699953642E-22, 0.0,
    -5.3143455936079533E-22, 0.0, -4.8638758239346971E-22, 0.0,
    -4.3021330714427634E-22, 0.0, -3.7069986106205997E-22, 0.0,
    -3.2411165937295763E-22, 0.0, -2.6670726717660928E-22, 0.0,
    -2.2255312440138803E-22, 0.0, -1.7604694431092617E-22, 0.0,
    -1.271851640695825E-22, 0.0, -8.410151184512858E-23, 0.0,
    -3.899483769772138E-23, 0.0, 0.0, 0.0 };

  static const real_T b_Kx[12] = { 7.689433750373557E-6, -0.00052983276131042512,
    -2.0904219222393924E-7, 0.033752729948155979, 6.2634718616689533E-6,
    -0.00049099767826871524, -1.9372005383038419E-7, 0.031251490572067622,
    5.0972730132940709E-6, -0.00045360951375905055, -1.7896878806682674E-7,
    0.028849976794694913 };

  static const real_T b_Ku1[3] = { 0.0056257498156150939, 0.0053016816587514205,
    0.0049805833709413629 };

  /* MATLAB Function: '<S1>/DataTypeConversion_reldist' incorporates:
   *  Inport: '<Root>/Relative distance'
   */
  Adaptive_DataTypeConversion_L0(Adaptive_U.d_rel, &rtb_y_ov);

  /* MATLAB Function: '<S1>/DataTypeConversion_vego' incorporates:
   *  Inport: '<Root>/Longitudinal velocity'
   */
  Adaptive_DataTypeConversion_L0(Adaptive_U.v_ego, &rtb_y_d);

  /* MATLAB Function: '<S1>/DataTypeConversion_L0' incorporates:
   *  Constant: '<S1>/Default spacing constant'
   */
  Adaptive_DataTypeConversion_L0(10.0, &rtb_y_j);

  /* MATLAB Function: '<S1>/DataTypeConversion_vset' incorporates:
   *  Inport: '<Root>/Set velocity'
   */
  Adaptive_DataTypeConversion_L0(Adaptive_U.set_velocity, &rtb_y);

  /* MATLAB Function: '<S1>/DataTypeConversion_vlead' incorporates:
   *  Inport: '<Root>/Longitudinal velocity'
   *  Inport: '<Root>/Relative velocity'
   *  Sum: '<S1>/Sum6'
   */
  Adaptive_DataTypeConversion_L0(Adaptive_U.v_ego + Adaptive_U.relative_velocity,
    &rtb_y_g);

  /* MATLAB Function: '<S1>/DataTypeConversion_amin' incorporates:
   *  Constant: '<S1>/Minimum longitudinal acceleration constant'
   */
  Adaptive_DataTypeConversion_L0(-3.0, &rtb_y_b);

  /* MATLAB Function: '<S1>/DataTypeConversion_amax' incorporates:
   *  Constant: '<S1>/Maximum longitudinal acceleration constant'
   */
  Adaptive_DataTypeConversion_L0(2.0, &rtb_y_ds);

  /* SignalConversion generated from: '<S34>/ SFunction ' incorporates:
   *  Constant: '<S1>/Default spacing constant'
   *  Constant: '<S1>/Minimum velocity constant'
   *  Inport: '<Root>/Longitudinal velocity'
   *  Inport: '<Root>/Time gap'
   *  MATLAB Function: '<S1>/DataTypeConversion_dmin'
   *  MATLAB Function: '<S33>/optimizer'
   *  Product: '<S1>/Product2'
   *  Sum: '<S1>/Sum1'
   */
  Adaptive_DataTypeConversion_L0(Adaptive_U.v_ego * Adaptive_U.Timegap + 10.0,
    &rtb_TmpSignalConversionAtSFun_e[0]);
  rtb_TmpSignalConversionAtSFun_e[1] = 0.0;

  /* SignalConversion generated from: '<S34>/ SFunction ' incorporates:
   *  Constant: '<S1>/Maximum velocity constant'
   *  Constant: '<S1>/Unconstrained'
   *  MATLAB Function: '<S33>/optimizer'
   */
  rtb_TmpSignalConversionAtSFu_mm[0] = 0;
  rtb_TmpSignalConversionAtSFu_mm[1] = 50;

  /* MATLAB Function: '<S33>/optimizer' incorporates:
   *  Memory: '<S13>/last_x'
   *  SignalConversion generated from: '<S34>/ SFunction '
   */
  memset(&vseq[0], 0, 62U * sizeof(real_T));
  for (i = 0; i < 31; i++) {
    vseq[(i << 1) + 1] = 1.0;
  }

  for (i = 0; i < 30; i++) {
    d_i = i << 1;
    rseq[d_i] = rtb_y_j * 0.02 - 0.76;
    rseq[d_i + 1] = rtb_y * 0.02 - 0.4;
  }

  for (i = 0; i < 31; i++) {
    vseq[i << 1] = Adaptive_RMDscale * rtb_y_g - Adaptive_voff;
  }

  rtb_y_j = vseq[0];
  rtb_y = vseq[1];
  xk[0] = Adaptive_DW.last_x_PreviousInput[0];
  xk[1] = Adaptive_DW.last_x_PreviousInput[1];
  xk[2] = Adaptive_DW.last_x_PreviousInput[2];
  xk[3] = Adaptive_DW.last_x_PreviousInput[3];

  /* SignalConversion generated from: '<S34>/ SFunction ' incorporates:
   *  MATLAB Function: '<S33>/optimizer'
   */
  ymax_incr[0] = rtb_y_ov * 0.02 - 0.76;
  ymax_incr[1] = rtb_y_d * 0.02 - 0.4;

  /* MATLAB Function: '<S33>/optimizer' incorporates:
   *  Memory: '<S13>/last_x'
   */
  rtb_y_ov = Adaptive_DW.last_x_PreviousInput[1];
  rtb_y_d = Adaptive_DW.last_x_PreviousInput[0];
  rtb_y_g = Adaptive_DW.last_x_PreviousInput[2];
  xk_0 = Adaptive_DW.last_x_PreviousInput[3];
  for (i = 0; i <= 0; i += 2) {
    /* MATLAB Function: '<S33>/optimizer' */
    tmp = _mm_loadu_pd(&ymax_incr[i]);
    _mm_storeu_pd(&y_innov[i], _mm_sub_pd(tmp, _mm_add_pd(_mm_add_pd(_mm_add_pd
      (_mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&f_a[i + 2]), _mm_set1_pd(rtb_y_ov)),
                  _mm_mul_pd(_mm_loadu_pd(&f_a[i]), _mm_set1_pd(rtb_y_d))),
       _mm_mul_pd(_mm_loadu_pd(&f_a[i + 4]), _mm_set1_pd(rtb_y_g))), _mm_mul_pd
      (_mm_loadu_pd(&f_a[i + 6]), _mm_set1_pd(xk_0))), _mm_set1_pd(0.0 * rtb_y_j
      + 0.0 * rtb_y))));
  }

  /* MATLAB Function: '<S33>/optimizer' */
  y_innov_0 = y_innov[1];
  y_innov_1 = y_innov[0];
  for (i = 0; i <= 2; i += 2) {
    /* MATLAB Function: '<S33>/optimizer' */
    tmp = _mm_loadu_pd(&xk[i]);
    _mm_storeu_pd(&rtb_xest[i], _mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_loadu_pd
      (&e_a[i + 4]), _mm_set1_pd(y_innov_0)), _mm_mul_pd(_mm_loadu_pd(&e_a[i]),
      _mm_set1_pd(y_innov_1))), tmp));
  }

  /* MATLAB Function: '<S33>/optimizer' incorporates:
   *  UnitDelay: '<S13>/last_mv'
   */
  ymax_incr_flag[0] = false;
  ymax_incr[0] = 0.0;
  ymin_incr_flag[0] = false;
  ymin_incr[0] = 0.0;
  ymax_incr_flag[1] = false;
  ymax_incr[1] = 0.0;
  ymin_incr_flag[1] = false;
  ymin_incr[1] = 0.0;
  umax_incr_flag = false;
  rtb_y_ov = 0.0;
  umin_incr_flag = false;
  rtb_y_d = 0.0;
  for (d_i = 0; d_i < 96; d_i++) {
    rtb_y_g = 0.0;
    for (i = 0; i < 62; i++) {
      rtb_y_g += b_Mv[96 * i + d_i] * vseq[i];
    }

    xk_0 = b_Mlim[d_i];
    rtb_y_g = -((((((b_Mx[d_i + 96] * rtb_xest[1] + b_Mx[d_i] * rtb_xest[0]) +
                    b_Mx[d_i + 192] * rtb_xest[2]) + b_Mx[d_i + 288] * rtb_xest
                   [3]) + xk_0) + b_Mu1[d_i] * Adaptive_DW.last_mv_DSTATE) +
                rtb_y_g);
    Bc[d_i] = rtb_y_g;
    b_Mrows = b_Mrows_0[d_i];
    if (b_Mrows <= 60) {
      i = (b_Mrows - (((b_Mrows - 1) / Adaptive_ny) << 1)) - 1;
      b_Del_Save_Flag0 = ymax_incr_flag[i];
      if (!ymax_incr_flag[i]) {
        xk_0 = -(0.02 * (real_T)rtb_TmpSignalConversionAtSFu_mm[i] - (-0.36 *
                  (real_T)i + 0.76)) - (-xk_0);
        b_Del_Save_Flag0 = true;
      } else {
        xk_0 = ymax_incr[i];
      }

      ymax_incr[i] = xk_0;
      ymax_incr_flag[i] = b_Del_Save_Flag0;
      Bc[d_i] = rtb_y_g + xk_0;
    } else if (b_Mrows <= 120) {
      i = (b_Mrows - (((b_Mrows - 61) >> 1) << 1)) - 61;
      b_Del_Save_Flag0 = ymin_incr_flag[i];
      if (!ymin_incr_flag[i]) {
        xk_0 = (0.02 * rtb_TmpSignalConversionAtSFun_e[i] - (-0.36 * (real_T)i +
                 0.76)) - (-xk_0);
        b_Del_Save_Flag0 = true;
      } else {
        xk_0 = ymin_incr[i];
      }

      ymin_incr[i] = xk_0;
      ymin_incr_flag[i] = b_Del_Save_Flag0;
      Bc[d_i] = rtb_y_g + xk_0;
    } else if (b_Mrows <= 150) {
      if (!umax_incr_flag) {
        rtb_y_ov = -(Adaptive_RMVscale * rtb_y_ds) - (-xk_0);
        umax_incr_flag = true;
      }

      Bc[d_i] = rtb_y_g + rtb_y_ov;
    } else {
      if (!umin_incr_flag) {
        rtb_y_d = Adaptive_RMVscale * rtb_y_b - (-xk_0);
        umin_incr_flag = true;
      }

      Bc[d_i] = rtb_y_g + rtb_y_d;
    }
  }

  f[0] = 0.0;
  f[1] = 0.0;
  f[2] = 0.0;
  f[3] = 0.0;
  for (d_i = 0; d_i < 3; d_i++) {
    rtb_y_b = 0.0;
    for (i = 0; i < 60; i++) {
      rtb_y_b += b_Kr[60 * d_i + i] * rseq[i];
    }

    rtb_y_ds = 0.0;
    for (i = 0; i < 62; i++) {
      rtb_y_ds += b_Kv[62 * d_i + i] * vseq[i];
    }

    i = d_i << 2;
    f[d_i] = (((((b_Kx[i + 1] * rtb_xest[1] + b_Kx[i] * rtb_xest[0]) + b_Kx[i +
                 2] * rtb_xest[2]) + b_Kx[i + 3] * rtb_xest[3]) + rtb_y_b) +
              b_Ku1[d_i] * Adaptive_DW.last_mv_DSTATE) + rtb_y_ds;
  }

  Adaptive_qpkwik(b_Linv, b_Hinv, f, b_Ac, Bc, Adaptive_DW.Memory_PreviousInput,
                  400, 1.0E-6, rtb_xest, a__1, &i);
  if ((i < 0) || (i == 0)) {
    rtb_xest[0] = 0.0;
  }

  Adaptive_DW.last_mv_DSTATE += rtb_xest[0];

  /* Outport: '<Root>/Longitudinal acceleration' incorporates:
   *  Gain: '<S13>/umin_scale1'
   *  MATLAB Function: '<S33>/optimizer'
   */
  Adaptive_Y.Longitudinalacceleration = 5.0 * Adaptive_DW.last_mv_DSTATE;

  /* MATLAB Function: '<S1>/DataTypeConversion_atrack' incorporates:
   *  Constant: '<S1>/External control signal constant'
   */
  Adaptive_DataTypeConversion_L0(0.0, &rtb_y_b);

  /* MATLAB Function: '<S33>/optimizer' */
  rtb_y_ov = xk[1];
  rtb_y_d = xk[0];
  rtb_y_g = xk[2];
  xk_0 = xk[3];
  y_innov_0 = y_innov[1];
  y_innov_1 = y_innov[0];
  for (i = 0; i <= 2; i += 2) {
    /* Update for Memory: '<S13>/last_x' incorporates:
     *  MATLAB Function: '<S33>/optimizer'
     */
    _mm_storeu_pd(&Adaptive_DW.last_x_PreviousInput[i], _mm_add_pd(_mm_add_pd
      (_mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(_mm_loadu_pd
      (&d_a[i + 4]), _mm_set1_pd(rtb_y_ov)), _mm_mul_pd(_mm_loadu_pd(&d_a[i]),
      _mm_set1_pd(rtb_y_d))), _mm_mul_pd(_mm_loadu_pd(&d_a[i + 8]), _mm_set1_pd
      (rtb_y_g))), _mm_mul_pd(_mm_loadu_pd(&d_a[i + 12]), _mm_set1_pd(xk_0))),
                  _mm_mul_pd(_mm_loadu_pd(&c_a[i]), _mm_set1_pd
      (Adaptive_DW.last_mv_DSTATE))), _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_a[i]),
      _mm_set1_pd(rtb_y_j)), _mm_set1_pd(0.0 * rtb_y))), _mm_add_pd(_mm_mul_pd
      (_mm_loadu_pd(&a[i + 4]), _mm_set1_pd(y_innov_0)), _mm_mul_pd(_mm_loadu_pd
      (&a[i]), _mm_set1_pd(y_innov_1)))));
  }
}

/* Model initialize function */
void Adaptive_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));
}

/* Model terminate function */
void Adaptive_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
