/*
 * Adaptive_sid.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Adaptive_sf".
 *
 * Model version              : 11.0
 * Simulink Coder version : 9.9 (R2023a) 19-Nov-2022
 * C source code generated on : Thu Aug  3 20:54:02 2023
 *
 * Target selection: rtwsfcn.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 *
 * SOURCES: Adaptive_sf.c
 */

/* statically allocated instance data for model: Adaptive */
{
  extern const ConstB_Adaptive_T Adaptive_Invariant;

  {
    /* Local SimStruct for the generated S-Function */
    static LocalS slS;
    LocalS *lS = &slS;
    ssSetUserData(rts, lS);

    /* block I/O */
    {
      static B_Adaptive_T sfcnB;
      void *b = (real_T *) &sfcnB;
      ssSetLocalBlockIO(rts, b);
      (void) memset(b, 0,
                    sizeof(B_Adaptive_T));
    }

    _ssSetConstBlockIO(rts, &Adaptive_Invariant);

    /* model checksums */
    ssSetChecksumVal(rts, 0, 1633255356U);
    ssSetChecksumVal(rts, 1, 2017059578U);
    ssSetChecksumVal(rts, 2, 1533458369U);
    ssSetChecksumVal(rts, 3, 2795744744U);
  }
}
