/** Four signal channels, each an array of intensities indexed by sample point. */
export interface ChannelSignals {
  A: number[];
  C: number[];
  G: number[];
  T: number[];
}

/**
 * Per-base chromatogram view: where each called base peaks within the signal
 * arrays, the four signal channels themselves, and optional per-base data.
 *
 * Generic on purpose — any concrete read model (e.g. `@fishka/seqio`'s parsed
 * chromatogram) is structurally compatible as long as it carries these fields.
 */
export interface Chromatogram {
  /** Per-base peak position as a sample-point index into the signal arrays. */
  positions: number[];
  /** A/C/G/T intensities by sample point; the four channels share one length. */
  signals: ChannelSignals;
  /** Per-base Phred-like confidence, one per called base. */
  confidences?: number[];
  /** Average peak spacing in sample points per base. */
  samplingRate?: number;
}
