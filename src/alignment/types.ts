/** Core generic types for the alignment utilities. Concrete project read/alignment
 * models (e.g. carrying chromatograms or extra metadata) are structurally
 * compatible as long as they include these fields. */

/** A half-open or inclusive span `[start, end]`; meaning is defined per function. */
export interface Region {
  start: number;
  end: number;
}

/** A named sequence (may contain gaps). */
export interface Sequence {
  name: string;
  sequence: string;
}

/** A sequence placed within a gapped alignment at a 0-based column offset. */
export interface AlignedSequence {
  /** Visual name of the sequence. */
  name: string;
  /** 0-based offset into the gapped reference where `sequence` starts. */
  offset: number;
  /** Sequence (with inner gaps). */
  sequence: string;
}

/** A read inside an alignment. Every read carries a stable `id`. */
export interface AlignedRead extends AlignedSequence {
  /** True if the read was reverse complemented to match reference orientation. */
  isReverseComplement: boolean;
  /** Stable per-read identifier, mandatory for every read inside an alignment. */
  id: string;
}

/** A reference plus its aligned reads. */
export interface Alignment {
  name: string;
  reference: AlignedSequence;
  reads: AlignedRead[];
}

/** One read aligned to a reference (both gapped). */
export interface PairwiseAlignment {
  /** Reference sequence with leading/inner/trailing gaps. */
  reference: Sequence;
  /** Read sequence with leading/inner/trailing gaps. */
  read: Sequence;
  /** True if the original read was reverse complemented. */
  isReverseComplement: boolean;
  /**
   * Stable read id, carried through so the resulting read keeps the same id
   * across (re)alignment. Optional here (a generic caller may omit it);
   * `decomposeToPairwiseAlignments` always populates it from the read.
   */
  id?: string;
}

/** A square alignment block: a gapped reference and equal-length gapped reads. */
export interface CombinedAlignment {
  /** A gapped reference sequence. */
  reference: Sequence;
  /** List of gapped read sequences. */
  reads: Sequence[];
}
