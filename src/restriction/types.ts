/** A restriction enzyme: its recognition sequence and where it cleaves each strand. */
export interface RestrictionEnzyme {
  /** Enzyme name (e.g. 'EcoRI'). */
  readonly name: string;
  /** Recognition sequence 5'->3' on the top strand, in IUPAC codes (e.g. 'GAATTC'). */
  readonly site: string;
  /**
   * Top-strand cut, as a bond index measured from the 5' end of {@link site}: `0` = before the first
   * base, `site.length` = after the last. EcoRI (G^AATTC) is `1`.
   */
  readonly cutTop: number;
  /**
   * Bottom-strand cut, in the SAME coordinate frame as {@link cutTop} (bases from the 5' end of the
   * top-strand recognition sequence). EcoRI is `5`. `cutBottom - cutTop` is the overhang: positive = a
   * 5' overhang, negative = a 3' overhang, `0` = a blunt cut.
   */
  readonly cutBottom: number;
}

/** One occurrence of an enzyme's recognition site on a sequence, resolved to absolute cut coordinates. */
export interface RestrictionSite {
  /** Enzyme name that produced this site. */
  readonly enzyme: string;
  /** The enzyme's recognition sequence (for reference / tooltip). */
  readonly site: string;
  /** 0-based start of the recognition match on the top strand. */
  readonly start: number;
  /** Recognition length in bases (`site.length`). */
  readonly length: number;
  /** Orientation of the recognition occurrence: `1` forward, `-1` on the reverse strand. Palindromic
   *  sites are always reported once, as `1`. */
  readonly strand: 1 | -1;
  /** Absolute top-strand cut, as a 0-based bond coordinate in `[0, sequenceLength]` (the primary tick
   *  position). For a circular molecule it is taken modulo the length. */
  readonly position: number;
  /** Absolute bottom-strand cut, same coordinate system as {@link position}. */
  readonly bottomCut: number;
  /** Overhang produced (`cutBottom - cutTop`): >0 a 5' overhang, <0 a 3' overhang, 0 blunt. */
  readonly overhang: number;
  /** True when the recognition site wraps across the origin of a circular molecule. */
  readonly crossesOrigin: boolean;
}
