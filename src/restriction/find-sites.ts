import { IUPAC_NUCLEOTIDES_LIST } from '../sequence/iupac.utils';
import { reverseComplement } from '../sequence/sequence.utils';
import { RestrictionEnzyme, RestrictionSite } from './types';

export interface FindRestrictionSitesOptions {
  /** Treat the molecule as circular: also find sites spanning the origin, and wrap cut coordinates. */
  readonly circular?: boolean;
}

/** Compiles a recognition string into a per-position set of the bases each IUPAC code admits. */
function compile(site: string): Set<string>[] {
  return [...site.toUpperCase()].map(ch => {
    const bases = IUPAC_NUCLEOTIDES_LIST[ch];
    if (!bases) throw new Error(`Invalid IUPAC code in recognition site: ${ch}`);
    return new Set(bases);
  });
}

/**
 * Subset match (same semantics as {@link findMotif}): a concrete read base matches iff the recognition
 * admits it; an ambiguous read base matches only if ALL the bases it stands for are admitted — so an `N`
 * run never spuriously reports a cut site.
 */
function baseMatches(base: string, allowed: Set<string>): boolean {
  const bases = IUPAC_NUCLEOTIDES_LIST[base?.toUpperCase()];
  if (!bases) return false; // gaps / unknown symbols match nothing
  for (const b of bases) if (!allowed.has(b)) return false;
  return true;
}

/** Whether a recognition site is its own reverse complement (so forward = reverse; report it once). */
function isPalindromic(site: string): boolean {
  const s = site.toUpperCase();
  return s === reverseComplement(s);
}

/**
 * Scan one strand's compiled recognition over the sequence, honouring circular wrap-around. Returns the
 * 0-based start index of every match; for a circular molecule a match may start near the end and read
 * through the origin (start stays in `[0, length)`).
 */
function scanStarts(sequence: string, motif: Set<string>[], circular: boolean): number[] {
  const length = sequence.length;
  const width = motif.length;
  if (width === 0 || width > length) return [];
  const lastStart = circular ? length - 1 : length - width;
  const starts: number[] = [];
  for (let i = 0; i <= lastStart; i++) {
    let hit = true;
    for (let j = 0; j < width; j++) {
      // Circular reads wrap; linear reads never index past `lastStart + width - 1 = length - 1`.
      const base = sequence[circular ? (i + j) % length : i + j];
      if (!baseMatches(base, motif[j])) {
        hit = false;
        break;
      }
    }
    if (hit) starts.push(i);
  }
  return starts;
}

/** Build a {@link RestrictionSite} from a match start, for a given orientation. */
function makeSite(
  enzyme: RestrictionEnzyme,
  start: number,
  strand: 1 | -1,
  length: number,
  circular: boolean,
): RestrictionSite {
  const L = enzyme.site.length;
  // Forward: cuts sit at `start + offset`. Reverse (recognition on the bottom strand): the enzyme's own
  // 5' end is at the high-coordinate end of the match, so its offsets mirror to `start + L - offset`,
  // and the two strands swap (the enzyme's top cut lands on our bottom strand and vice-versa).
  const rawTop = strand === 1 ? start + enzyme.cutTop : start + L - enzyme.cutBottom;
  const rawBottom = strand === 1 ? start + enzyme.cutBottom : start + L - enzyme.cutTop;
  const wrap = (n: number): number => (circular ? ((n % length) + length) % length : n);
  return {
    enzyme: enzyme.name,
    site: enzyme.site,
    start,
    length: L,
    strand,
    position: wrap(rawTop),
    bottomCut: wrap(rawBottom),
    overhang: enzyme.cutBottom - enzyme.cutTop,
    crossesOrigin: circular && start + L > length,
  };
}

/**
 * Find every cut site of the given enzymes in `sequence`. Both strands are searched for a
 * non-palindromic recognition sequence; a palindromic one is reported once (as `strand: 1`) since its
 * forward and reverse occurrences coincide. Cut coordinates are absolute top/bottom-strand bond
 * positions ({@link RestrictionSite.position} / {@link RestrictionSite.bottomCut}); on a circular
 * molecule ({@link FindRestrictionSitesOptions.circular}) sites spanning the origin are found and all
 * coordinates wrap into `[0, length)`.
 *
 * Results are sorted by top-strand cut position, then enzyme name, then strand. Throws when an enzyme's
 * recognition sequence contains a non-IUPAC character.
 *
 * Assumes a palindromic recognition sequence is cut symmetrically (`cutBottom === site.length -
 * cutTop`), which holds for every real Type II enzyme with a self-complementary site. A palindromic
 * site described with asymmetric cuts is not a valid enzyme and would be reported with only its
 * forward-orientation cut.
 */
export function findRestrictionSites(
  sequence: string,
  enzymes: readonly RestrictionEnzyme[],
  options: FindRestrictionSitesOptions = {},
): RestrictionSite[] {
  const circular = options.circular === true;
  const length = sequence.length;
  const sites: RestrictionSite[] = [];

  for (const enzyme of enzymes) {
    const forward = compile(enzyme.site);
    for (const start of scanStarts(sequence, forward, circular)) {
      sites.push(makeSite(enzyme, start, 1, length, circular));
    }
    if (isPalindromic(enzyme.site)) continue;
    const reverse = compile(reverseComplement(enzyme.site));
    for (const start of scanStarts(sequence, reverse, circular)) {
      sites.push(makeSite(enzyme, start, -1, length, circular));
    }
  }

  return sites.sort(
    (a, b) =>
      a.position - b.position || (a.enzyme < b.enzyme ? -1 : a.enzyme > b.enzyme ? 1 : 0) || a.strand - b.strand,
  );
}
