import { applyAlignmentGapsToChromatogram, Chromatogram, reverseComplementChromatogram } from '../../src/chromatogram';

/** Build a chromatogram whose every channel is sequential 0..signalLength-1, like the source tests. */
function createChromatogram(positions: number[], signalLength: number): Chromatogram {
  const channel = Array.from({ length: signalLength }, (_, i) => i);
  return {
    positions,
    signals: { A: [...channel], C: [...channel], G: [...channel], T: [...channel] },
    confidences: positions.map((_, i) => i % 100),
    samplingRate: 1,
  };
}

describe('reverseComplementChromatogram', () => {
  // Distinct values per channel so a wrong swap or missed reverse is caught exactly.
  // Asymmetric positions (not palindromic) so the L-1-p mirror is observable.
  const chromatogram: Chromatogram = {
    positions: [0, 1, 3],
    signals: {
      A: [1, 2, 3, 4, 5],
      C: [10, 20, 30, 40, 50],
      G: [11, 21, 31, 41, 51],
      T: [100, 200, 300, 400, 500],
    },
    confidences: [9, 8, 7],
    samplingRate: 2,
  };

  it('swaps A<->T and C<->G and reverses each channel', () => {
    const result = reverseComplementChromatogram(chromatogram);
    expect(result.signals.A).toEqual([500, 400, 300, 200, 100]); // reverse of old T
    expect(result.signals.T).toEqual([5, 4, 3, 2, 1]); // reverse of old A
    expect(result.signals.C).toEqual([51, 41, 31, 21, 11]); // reverse of old G
    expect(result.signals.G).toEqual([50, 40, 30, 20, 10]); // reverse of old C
  });

  it('mirrors positions against the signal length and re-sorts ascending', () => {
    // L = 5; [0,1,3] -> map(p => 4 - p) = [4,3,1] -> reverse = [1,3,4]
    expect(reverseComplementChromatogram(chromatogram).positions).toEqual([1, 3, 4]);
  });

  it('reverses confidences and keeps samplingRate', () => {
    const result = reverseComplementChromatogram(chromatogram);
    expect(result.confidences).toEqual([7, 8, 9]);
    expect(result.samplingRate).toBe(2);
  });

  it('leaves confidences undefined when there are none', () => {
    const result = reverseComplementChromatogram({
      positions: [0],
      signals: { A: [1], C: [2], G: [3], T: [4] },
    });
    expect(result.confidences).toBeUndefined();
  });

  it('is its own inverse (applying it twice returns the original)', () => {
    expect(reverseComplementChromatogram(reverseComplementChromatogram(chromatogram))).toEqual(chromatogram);
  });

  it('does not mutate the input', () => {
    const snapshot = JSON.parse(JSON.stringify(chromatogram));
    reverseComplementChromatogram(chromatogram);
    expect(chromatogram).toEqual(snapshot);
  });
});

describe('applyAlignmentGapsToChromatogram', () => {
  it('inserts leading gaps', () => {
    const result = applyAlignmentGapsToChromatogram(createChromatogram([5, 15, 25, 35], 40), '--ATGC');
    expect(result.positions).toEqual([0, 1, 7, 17, 27, 37]);
    expect(result.confidences).toEqual([0, 0, 0, 1, 2, 3]);
    expect(result.signals.A).toHaveLength(42);
    expect(result.signals.A[0]).toBe(0); // inserted gap sample
    expect(result.signals.A[2]).toBe(0); // original A[0], shifted by 2
    expect(result.signals.A[7]).toBe(5); // original A[5]
    expect(result.signals.A[41]).toBe(39); // last original sample
  });

  it('inserts internal gaps', () => {
    const result = applyAlignmentGapsToChromatogram(createChromatogram([5, 15, 25], 30), 'A--TG');
    expect(result.positions).toEqual([5, 15, 16, 17, 27]);
    expect(result.confidences).toEqual([0, 0, 0, 1, 2]);
    expect(result.signals.A).toHaveLength(32);
    expect(result.signals.A[15]).toBe(0); // first inserted gap sample
    expect(result.signals.A[17]).toBe(15); // original A[15], shifted by 2
  });

  it('inserts leading and internal gaps together', () => {
    const result = applyAlignmentGapsToChromatogram(createChromatogram([5, 15, 25], 30), '-A-T-G');
    expect(result.positions).toEqual([0, 6, 16, 17, 27, 28]);
    expect(result.confidences).toEqual([0, 0, 0, 1, 0, 2]);
    expect(result.signals.A).toHaveLength(33);
  });

  it('inserts consecutive internal gaps', () => {
    const result = applyAlignmentGapsToChromatogram(createChromatogram([5, 15], 20), 'A---T');
    expect(result.positions).toEqual([5, 15, 16, 17, 18]);
    expect(result.confidences).toEqual([0, 0, 0, 0, 1]);
    expect(result.signals.A).toHaveLength(23);
  });

  it('inserts trailing gaps', () => {
    const result = applyAlignmentGapsToChromatogram(createChromatogram([5, 15], 20), 'AT--');
    expect(result.positions).toEqual([5, 15, 20, 21]);
    expect(result.confidences).toEqual([0, 1, 0, 0]);
    expect(result.signals.A).toHaveLength(22);
  });

  it('returns the same chromatogram unchanged when there are no gaps', () => {
    const chromatogram = createChromatogram([5, 15, 25, 35], 40);
    expect(applyAlignmentGapsToChromatogram(chromatogram, 'ATGC')).toBe(chromatogram);
  });

  it('throws when the chromatogram length does not match the clean sequence', () => {
    expect(() => applyAlignmentGapsToChromatogram(createChromatogram([5, 15], 20), 'A-TGC')).toThrow();
  });
});
