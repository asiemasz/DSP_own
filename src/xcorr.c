/**
 * It's int8_t version of CMSIS DSP's correlate function for float32 type.
 */

#include "xcorr.h"
#include <assert.h>

void xcorrelate_int8(int8_t *data1, uint16_t data1Length, int8_t *data2,
                     uint16_t data2Length, int8_t *corr, uint16_t corrLength) {

  int8_t *pIn1;
  int8_t *pIn2;
  int8_t *pOut = corr;
  int8_t *px;
  int8_t *py;
  int8_t *pSrc1;
  int8_t sum, acc0, acc1, acc2, acc3;
  int32_t x0, x1, x2, x3, c0;
  uint32_t j, k = 0U, count, blkCnt, outBlockSize, blockSize1, blockSize2,
              blockSize3; /* loop counters */
  int32_t inc = 1;        /* Destination address modifier */

  if (data1Length >= data2Length) {
    pIn1 = data1;
    pIn2 = data2;

    outBlockSize = (2U * data1Length) - 1U;
    assert(corrLength == outBlockSize);

    // Zero padding for data2
    j = outBlockSize - (data1Length + (data2Length - 1U));

    pOut += j;
  } else {
    pIn1 = data2;
    pIn2 = data1;

    j = data2Length;
    data2Length = data1Length;
    data1Length = j;

    pOut = corr + ((data1Length + data2Length) - 2U);

    inc = -1;
  }

  blockSize1 = data2Length - 1U;
  blockSize2 = data1Length - (data2Length - 1U);
  blockSize3 = blockSize1;

  count = 1U;

  px = pIn1;

  pSrc1 = pIn2 + (data2Length - 1U);
  py = pSrc1;

  /* ------------------------
   * Stage1 process
   * ----------------------*/

  /* The first stage starts here */
  while (blockSize1 > 0U) {
    /* Accumulator is made zero for every iteration */
    sum = 0;

    /* Apply loop unrolling and compute 4 MACs simultaneously. */
    k = count >> 2U;

    /* First part of the processing with loop unrolling.  Compute 4 MACs at a
     *time.
     ** a second loop below computes MACs for the remaining 1 to 3 samples. */
    while (k > 0U) {
      /* x[0] * y[data2Length - 4] */
      sum += *px++ * *py++;
      /* x[1] * y[data2Length - 3] */
      sum += *px++ * *py++;
      /* x[2] * y[data2Length - 2] */
      sum += *px++ * *py++;
      /* x[3] * y[data2Length - 1] */
      sum += *px++ * *py++;

      /* Decrement the loop counter */
      k--;
    }

    /* If the count is not a multiple of 4, compute any remaining MACs here.
     ** No loop unrolling is used. */
    k = count % 0x4U;

    while (k > 0U) {
      /* Perform the multiply-accumulate */
      /* x[0] * y[data2Length - 1] */
      sum += *px++ * *py++;

      /* Decrement the loop counter */
      k--;
    }

    /* Store the result in the accumulator in the destination buffer. */
    *pOut = sum;
    /* Destination pointer is updated according to the address modifier, inc */
    pOut += inc;

    /* Update the inputA and inputB pointers for next MAC calculation */
    py = pSrc1 - count;
    px = pIn1;

    /* Increment the MAC count */
    count++;

    /* Decrement the loop counter */
    blockSize1--;
  }

  /* --------------------------
   * Initializations of stage2
   * ------------------------*/

  /* sum = x[0] * y[0] + x[1] * y[1] +...+ x[data2Length-1] * y[data2Length-1]
   * sum = x[1] * y[0] + x[2] * y[1] +...+ x[data2Length] * y[data2Length-1]
   * ....
   * sum = x[data1Length-data2Length-2] * y[0] + x[data1Length-data2Length-1] *
   * y[1] +...+ x[data1Length-1] * y[data2Length-1]
   */

  /* Working pointer of inputA */
  px = pIn1;

  /* Working pointer of inputB */
  py = pIn2;

  /* count is index by which the pointer pIn1 to be incremented */
  count = 0U;

  /* -------------------
   * Stage2 process
   * ------------------*/

  /* Stage2 depends on data2Length as in this stage data2Length number of MACS
   * are performed. So, to loop unroll over blockSize2, data2Length should be
   * greater than or equal to 4, to loop unroll the data2Length loop */
  if (data2Length >= 4U) {
    /* Loop unroll over blockSize2, by 4 */
    blkCnt = blockSize2 >> 2U;

    while (blkCnt > 0U) {
      /* Set all accumulators to zero */
      acc0 = 0;
      acc1 = 0;
      acc2 = 0;
      acc3 = 0;

      /* read x[0], x[1], x[2] samples */
      x0 = *(px++);
      x1 = *(px++);
      x2 = *(px++);

      /* Apply loop unrolling and compute 4 MACs simultaneously. */
      k = data2Length >> 2U;

      /* First part of the processing with loop unrolling.  Compute 4 MACs at a
       *time.
       ** a second loop below computes MACs for the remaining 1 to 3 samples. */
      do {
        /* Read y[0] sample */
        c0 = *(py++);

        /* Read x[3] sample */
        x3 = *(px++);

        /* Perform the multiply-accumulate */
        /* acc0 +=  x[0] * y[0] */
        acc0 += x0 * c0;
        /* acc1 +=  x[1] * y[0] */
        acc1 += x1 * c0;
        /* acc2 +=  x[2] * y[0] */
        acc2 += x2 * c0;
        /* acc3 +=  x[3] * y[0] */
        acc3 += x3 * c0;

        /* Read y[1] sample */
        c0 = *(py++);

        /* Read x[4] sample */
        x0 = *(px++);

        /* Perform the multiply-accumulate */
        /* acc0 +=  x[1] * y[1] */
        acc0 += x1 * c0;
        /* acc1 +=  x[2] * y[1] */
        acc1 += x2 * c0;
        /* acc2 +=  x[3] * y[1] */
        acc2 += x3 * c0;
        /* acc3 +=  x[4] * y[1] */
        acc3 += x0 * c0;

        /* Read y[2] sample */
        c0 = *(py++);

        /* Read x[5] sample */
        x1 = *(px++);

        /* Perform the multiply-accumulates */
        /* acc0 +=  x[2] * y[2] */
        acc0 += x2 * c0;
        /* acc1 +=  x[3] * y[2] */
        acc1 += x3 * c0;
        /* acc2 +=  x[4] * y[2] */
        acc2 += x0 * c0;
        /* acc3 +=  x[5] * y[2] */
        acc3 += x1 * c0;

        /* Read y[3] sample */
        c0 = *(py++);

        /* Read x[6] sample */
        x2 = *(px++);

        /* Perform the multiply-accumulates */
        /* acc0 +=  x[3] * y[3] */
        acc0 += x3 * c0;
        /* acc1 +=  x[4] * y[3] */
        acc1 += x0 * c0;
        /* acc2 +=  x[5] * y[3] */
        acc2 += x1 * c0;
        /* acc3 +=  x[6] * y[3] */
        acc3 += x2 * c0;

      } while (--k);

      /* If the data2Length is not a multiple of 4, compute any remaining MACs
       *here.
       ** No loop unrolling is used. */
      k = data2Length % 0x4U;

      while (k > 0U) {
        /* Read y[4] sample */
        c0 = *(py++);

        /* Read x[7] sample */
        x3 = *(px++);

        /* Perform the multiply-accumulates */
        /* acc0 +=  x[4] * y[4] */
        acc0 += x0 * c0;
        /* acc1 +=  x[5] * y[4] */
        acc1 += x1 * c0;
        /* acc2 +=  x[6] * y[4] */
        acc2 += x2 * c0;
        /* acc3 +=  x[7] * y[4] */
        acc3 += x3 * c0;

        /* Reuse the present samples for the next MAC */
        x0 = x1;
        x1 = x2;
        x2 = x3;

        /* Decrement the loop counter */
        k--;
      }

      /* Store the result in the accumulator in the destination buffer. */
      *pOut = acc0;
      /* Destination pointer is updated according to the address modifier, inc
       */
      pOut += inc;

      *pOut = acc1;
      pOut += inc;

      *pOut = acc2;
      pOut += inc;

      *pOut = acc3;
      pOut += inc;

      /* Increment the pointer pIn1 index, count by 4 */
      count += 4U;

      /* Update the inputA and inputB pointers for next MAC calculation */
      px = pIn1 + count;
      py = pIn2;

      /* Decrement the loop counter */
      blkCnt--;
    }

    /* If the blockSize2 is not a multiple of 4, compute any remaining output
     *samples here.
     ** No loop unrolling is used. */
    blkCnt = blockSize2 % 0x4U;

    while (blkCnt > 0U) {
      /* Accumulator is made zero for every iteration */
      sum = 0;

      /* Apply loop unrolling and compute 4 MACs simultaneously. */
      k = data2Length >> 2U;

      /* First part of the processing with loop unrolling.  Compute 4 MACs at a
       *time.
       ** a second loop below computes MACs for the remaining 1 to 3 samples. */
      while (k > 0U) {
        /* Perform the multiply-accumulates */
        sum += *px++ * *py++;
        sum += *px++ * *py++;
        sum += *px++ * *py++;
        sum += *px++ * *py++;

        /* Decrement the loop counter */
        k--;
      }

      /* If the data2Length is not a multiple of 4, compute any remaining MACs
       *here.
       ** No loop unrolling is used. */
      k = data2Length % 0x4U;

      while (k > 0U) {
        /* Perform the multiply-accumulate */
        sum += *px++ * *py++;

        /* Decrement the loop counter */
        k--;
      }

      /* Store the result in the accumulator in the destination buffer. */
      *pOut = sum;
      /* Destination pointer is updated according to the address modifier, inc
       */
      pOut += inc;

      /* Increment the pointer pIn1 index, count by 1 */
      count++;

      /* Update the inputA and inputB pointers for next MAC calculation */
      px = pIn1 + count;
      py = pIn2;

      /* Decrement the loop counter */
      blkCnt--;
    }
  } else {
    /* If the data2Length is not a multiple of 4,
     * the blockSize2 loop cannot be unrolled by 4 */
    blkCnt = blockSize2;

    while (blkCnt > 0U) {
      /* Accumulator is made zero for every iteration */
      sum = 0;

      /* Loop over data2Length */
      k = data2Length;

      while (k > 0U) {
        /* Perform the multiply-accumulate */
        sum += *px++ * *py++;

        /* Decrement the loop counter */
        k--;
      }

      /* Store the result in the accumulator in the destination buffer. */
      *pOut = sum;
      /* Destination pointer is updated according to the address modifier, inc
       */
      pOut += inc;

      /* Increment the pointer pIn1 index, count by 1 */
      count++;

      /* Update the inputA and inputB pointers for next MAC calculation */
      px = pIn1 + count;
      py = pIn2;

      /* Decrement the loop counter */
      blkCnt--;
    }
  }

  /* --------------------------
   * Initializations of stage3
   * -------------------------*/

  /* sum += x[data1Length-data2Length+1] * y[0] + x[data1Length-data2Length+2] *
   * y[1] +...+ x[data1Length-1] * y[data2Length-1] sum +=
   * x[data1Length-data2Length+2] * y[0] + x[data1Length-data2Length+3] * y[1]
   * +...+ x[data1Length-1] * y[data2Length-1]
   * ....
   * sum +=  x[data1Length-2] * y[0] + x[data1Length-1] * y[1]
   * sum +=  x[data1Length-1] * y[0]
   */

  /* In this stage the MAC operations are decreased by 1 for every iteration.
     The count variable holds the number of MAC operations performed */
  count = data2Length - 1U;

  /* Working pointer of inputA */
  pSrc1 = pIn1 + (data1Length - (data2Length - 1U));
  px = pSrc1;

  /* Working pointer of inputB */
  py = pIn2;

  /* -------------------
   * Stage3 process
   * ------------------*/

  while (blockSize3 > 0U) {
    /* Accumulator is made zero for every iteration */
    sum = 0;

    /* Apply loop unrolling and compute 4 MACs simultaneously. */
    k = count >> 2U;

    /* First part of the processing with loop unrolling.  Compute 4 MACs at a
     *time.
     ** a second loop below computes MACs for the remaining 1 to 3 samples. */
    while (k > 0U) {
      /* Perform the multiply-accumulates */
      /* sum += x[data1Length - data2Length + 4] * y[3] */
      sum += *px++ * *py++;
      /* sum += x[data1Length - data2Length + 3] * y[2] */
      sum += *px++ * *py++;
      /* sum += x[data1Length - data2Length + 2] * y[1] */
      sum += *px++ * *py++;
      /* sum += x[data1Length - data2Length + 1] * y[0] */
      sum += *px++ * *py++;

      /* Decrement the loop counter */
      k--;
    }

    /* If the count is not a multiple of 4, compute any remaining MACs here.
     ** No loop unrolling is used. */
    k = count % 0x4U;

    while (k > 0U) {
      /* Perform the multiply-accumulates */
      sum += *px++ * *py++;

      /* Decrement the loop counter */
      k--;
    }

    /* Store the result in the accumulator in the destination buffer. */
    *pOut = sum;
    /* Destination pointer is updated according to the address modifier, inc */
    pOut += inc;

    /* Update the inputA and inputB pointers for next MAC calculation */
    px = ++pSrc1;
    py = pIn2;

    /* Decrement the MAC count */
    count--;

    /* Decrement the loop counter */
    blockSize3--;
  }
}