import { MottTrimOptions, SlidingWindowTrimOptions, TrimRange } from './types';

/** Kept-region default cutoff / threshold in Phred quality units (≈ error rate 0.01). */
const DEFAULT_QUALITY = 20;
/** Default sliding window width in bases. */
const DEFAULT_WINDOW_SIZE = 10;

const EMPTY_RANGE: TrimRange = { start: 0, end: 0 };

/**
 * Convert a Phred quality score `Q` to its error probability `10^(-Q/10)`.
 * Q20 → 0.01, Q30 → 0.001. Used to weight bases on the same scale phred does.
 */
export function phredToErrorProbability(quality: number): number {
  return Math.pow(10, -quality / 10);
}

/**
 * Trim a read to its highest-quality core with the **modified Mott algorithm**
 * (the one phred/Geneious use). Each base scores `errorProb(cutoff) −
 * errorProb(Q)`: bases cleaner than the cutoff score positive, dirtier ones
 * negative. The maximum-sum contiguous run of those scores is the kept region,
 * so noisy leading/trailing bases are dropped while a good core — even one that
 * dips briefly in the middle — is preserved.
 *
 * `qualities` is the per-base Phred confidence array (an .ab1 `PCON` tag).
 * Returns a half-open {@link TrimRange}; `{ start: 0, end: 0 }` when no run
 * scores positive (whole read below the cutoff) or the best run is shorter than
 * `minLength`.
 */
export function mottTrim(qualities: number[], options?: MottTrimOptions): TrimRange {
  const cutoff = options?.cutoff ?? DEFAULT_QUALITY;
  const minLength = options?.minLength ?? 0;
  const cutoffError = phredToErrorProbability(cutoff);

  // Kadane's max-subarray over per-base scores. bestSum starts at 0, so a region
  // is only chosen when its summed score is strictly positive.
  let bestStart = 0;
  let bestEnd = 0;
  let bestSum = 0;
  let runStart = 0;
  let runSum = 0;
  for (let i = 0; i < qualities.length; i++) {
    const score = cutoffError - phredToErrorProbability(qualities[i]);
    if (runSum <= 0) {
      runStart = i;
      runSum = score;
    } else {
      runSum += score;
    }
    if (runSum > bestSum) {
      bestSum = runSum;
      bestStart = runStart;
      bestEnd = i + 1;
    }
  }

  if (bestEnd - bestStart < minLength) {
    return EMPTY_RANGE;
  }
  return { start: bestStart, end: bestEnd };
}

/**
 * Trim a read by **sliding-window mean quality** (the Trimmomatic-style scan).
 * A window of `windowSize` bases is "good" when its mean Phred quality reaches
 * `threshold`. The kept region spans from the first good window to the last,
 * then each edge is pulled inward past any residual sub-threshold bases the
 * window dragged in, giving a tight boundary; interior dips a good window still
 * averages over are left intact. The window is clamped to the read length so
 * reads shorter than `windowSize` are still evaluated as a single window.
 *
 * `qualities` is the per-base Phred confidence array (an .ab1 `PCON` tag).
 * Returns a half-open {@link TrimRange}; `{ start: 0, end: 0 }` when no window
 * meets the threshold or the kept region is shorter than `minLength`.
 */
export function slidingWindowTrim(qualities: number[], options?: SlidingWindowTrimOptions): TrimRange {
  const threshold = options?.threshold ?? DEFAULT_QUALITY;
  const minLength = options?.minLength ?? 0;
  const requestedWindow = options?.windowSize ?? DEFAULT_WINDOW_SIZE;
  const n = qualities.length;
  if (n === 0) {
    return EMPTY_RANGE;
  }
  const window = Math.min(requestedWindow, n);

  let firstLeft = -1;
  let lastRight = -1;
  let sum = 0;
  for (let i = 0; i < window; i++) {
    sum += qualities[i];
  }
  for (let left = 0; left + window <= n; left++) {
    if (left > 0) {
      sum += qualities[left + window - 1] - qualities[left - 1];
    }
    if (sum / window >= threshold) {
      if (firstLeft === -1) {
        firstLeft = left;
      }
      lastRight = left + window;
    }
  }

  if (firstLeft === -1) {
    return EMPTY_RANGE;
  }

  // Pull each edge inward past sub-threshold bases the boundary windows averaged over.
  let start = firstLeft;
  let end = lastRight;
  while (start < end && qualities[start] < threshold) {
    start++;
  }
  while (end > start && qualities[end - 1] < threshold) {
    end--;
  }

  if (end - start < minLength) {
    return EMPTY_RANGE;
  }
  return { start, end };
}
