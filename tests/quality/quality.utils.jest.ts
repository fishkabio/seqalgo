import { mottTrim, phredToErrorProbability, slidingWindowTrim } from '../../src/quality';

/** A run of `n` bases all at quality `q`. */
function run(q: number, n: number): number[] {
  return new Array<number>(n).fill(q);
}

describe('phredToErrorProbability', () => {
  it('maps standard Phred scores to error probabilities', () => {
    expect(phredToErrorProbability(0)).toBeCloseTo(1);
    expect(phredToErrorProbability(10)).toBeCloseTo(0.1);
    expect(phredToErrorProbability(20)).toBeCloseTo(0.01);
    expect(phredToErrorProbability(30)).toBeCloseTo(0.001);
  });
});

describe('mottTrim', () => {
  it('drops low-quality flanks and keeps the good core', () => {
    // 5 dirty (Q5) + 20 clean (Q40) + 5 dirty (Q5)
    const qualities = [...run(5, 5), ...run(40, 20), ...run(5, 5)];
    expect(mottTrim(qualities)).toEqual({ start: 5, end: 25 });
  });

  it('keeps the whole read when every base is above the cutoff', () => {
    const qualities = run(40, 30);
    expect(mottTrim(qualities)).toEqual({ start: 0, end: 30 });
  });

  it('returns an empty range when the whole read is below the cutoff', () => {
    expect(mottTrim(run(5, 30))).toEqual({ start: 0, end: 0 });
  });

  it('returns an empty range for empty input', () => {
    expect(mottTrim([])).toEqual({ start: 0, end: 0 });
  });

  it('spans a shallow interior dip rather than splitting the read', () => {
    // A short, shallow dip (Q13, just below the cutoff) inside a long clean read
    // costs less than the flanking gain, so the kept region bridges it.
    const qualities = [...run(40, 15), ...run(13, 1), ...run(40, 15)];
    expect(mottTrim(qualities)).toEqual({ start: 0, end: 31 });
  });

  it('splits at a deep interior dip too costly to span', () => {
    // A deep dip (Q2) costs more than the flanking gain, so only the first
    // high-quality run is kept (the earliest maximal run wins ties).
    const qualities = [...run(40, 15), ...run(2, 2), ...run(40, 15)];
    expect(mottTrim(qualities)).toEqual({ start: 0, end: 15 });
  });

  it('respects a custom cutoff (higher cutoff trims mid-quality flanks)', () => {
    // Q15 flanks survive the default (20 → they score negative but... ) — check both cutoffs.
    const qualities = [...run(15, 5), ...run(40, 20), ...run(15, 5)];
    // Default cutoff 20: Q15 is below cutoff, so flanks are trimmed.
    expect(mottTrim(qualities)).toEqual({ start: 5, end: 25 });
    // Cutoff 10: Q15 is now above cutoff and kept.
    expect(mottTrim(qualities, { cutoff: 10 })).toEqual({ start: 0, end: 30 });
  });

  it('rejects a kept region shorter than minLength', () => {
    const qualities = [...run(5, 10), ...run(40, 3), ...run(5, 10)];
    expect(mottTrim(qualities, { minLength: 5 })).toEqual({ start: 0, end: 0 });
    expect(mottTrim(qualities, { minLength: 3 })).toEqual({ start: 10, end: 13 });
  });

  it('does not mutate the input', () => {
    const qualities = [...run(5, 3), ...run(40, 5)];
    const copy = [...qualities];
    mottTrim(qualities);
    expect(qualities).toEqual(copy);
  });
});

describe('slidingWindowTrim', () => {
  it('drops low-quality flanks and keeps the good core', () => {
    const qualities = [...run(5, 8), ...run(40, 20), ...run(5, 8)];
    expect(slidingWindowTrim(qualities, { windowSize: 4, threshold: 20 })).toEqual({ start: 8, end: 28 });
  });

  it('keeps the whole read when every window is above threshold', () => {
    const qualities = run(40, 30);
    expect(slidingWindowTrim(qualities, { windowSize: 10, threshold: 20 })).toEqual({ start: 0, end: 30 });
  });

  it('returns an empty range when no window meets the threshold', () => {
    expect(slidingWindowTrim(run(5, 30), { windowSize: 10, threshold: 20 })).toEqual({ start: 0, end: 0 });
  });

  it('returns an empty range for empty input', () => {
    expect(slidingWindowTrim([])).toEqual({ start: 0, end: 0 });
  });

  it('clamps the window to the read length for short reads', () => {
    // 3 bases, window requested 10 → evaluated as one window of 3; mean 40 ≥ 20.
    expect(slidingWindowTrim(run(40, 3), { windowSize: 10, threshold: 20 })).toEqual({ start: 0, end: 3 });
  });

  it('keeps an interior dip that a good window still averages over', () => {
    // A 2-base dip inside a window of 4 with Q40 neighbours: window mean stays ≥ 20.
    const qualities = [...run(40, 6), ...run(2, 2), ...run(40, 6)];
    expect(slidingWindowTrim(qualities, { windowSize: 4, threshold: 20 })).toEqual({ start: 0, end: 14 });
  });

  it('rejects a kept region shorter than minLength', () => {
    const qualities = [...run(5, 10), ...run(40, 4), ...run(5, 10)];
    expect(slidingWindowTrim(qualities, { windowSize: 4, threshold: 20, minLength: 5 })).toEqual({ start: 0, end: 0 });
    expect(slidingWindowTrim(qualities, { windowSize: 4, threshold: 20, minLength: 4 })).toEqual({ start: 10, end: 14 });
  });

  it('does not mutate the input', () => {
    const qualities = [...run(5, 4), ...run(40, 6)];
    const copy = [...qualities];
    slidingWindowTrim(qualities, { windowSize: 4 });
    expect(qualities).toEqual(copy);
  });
});
