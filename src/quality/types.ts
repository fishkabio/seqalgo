/**
 * A half-open kept-region over a per-base quality array: bases `[start, end)`
 * are retained, everything outside is trimmed. `start === end` means "trim
 * everything" (no base survived the quality criterion).
 */
export interface TrimRange {
  /** Index of the first kept base (inclusive). */
  start: number;
  /** Index one past the last kept base (exclusive). */
  end: number;
}

/** Options for {@link mottTrim}. */
export interface MottTrimOptions {
  /**
   * Phred quality cutoff. A base contributes positively to the retained region
   * when its quality exceeds this and negatively when it falls below, so the
   * boundaries land where good signal starts and stops. Default 20.
   */
  cutoff?: number;
  /**
   * Minimum length of the kept region. A shorter result is rejected as an empty
   * range (`{ start: 0, end: 0 }`). Default 0 (no minimum).
   */
  minLength?: number;
}

/** Options for {@link slidingWindowTrim}. */
export interface SlidingWindowTrimOptions {
  /** Window width in bases. Default 10. Clamped to the read length for short reads. */
  windowSize?: number;
  /** Mean Phred quality a window must reach to count as "good". Default 20. */
  threshold?: number;
  /**
   * Minimum length of the kept region. A shorter result is rejected as an empty
   * range (`{ start: 0, end: 0 }`). Default 0 (no minimum).
   */
  minLength?: number;
}
