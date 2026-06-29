import { Chromatogram, reverseComplementChromatogram } from '../../src/chromatogram';

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
