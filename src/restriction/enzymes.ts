import { RestrictionEnzyme } from './types';

/**
 * A curated set of ~25 common cloning restriction enzymes (all Type II, palindromic recognition).
 * Each entry carries the recognition sequence 5'->3' and the two cut offsets (see
 * {@link RestrictionEnzyme}). Values follow the standard REBASE cut notation, e.g. EcoRI G^AATT_C →
 * cutTop 1, cutBottom 5 (a 4-nt 5' overhang).
 *
 * This is a deliberately small, license-clean starter set — not a full REBASE mirror. Callers can pass
 * their own {@link RestrictionEnzyme} list to {@link findRestrictionSites} for anything beyond it.
 */
export const COMMON_RESTRICTION_ENZYMES: readonly RestrictionEnzyme[] = [
  // 5' overhangs.
  { name: 'EcoRI', site: 'GAATTC', cutTop: 1, cutBottom: 5 },
  { name: 'BamHI', site: 'GGATCC', cutTop: 1, cutBottom: 5 },
  { name: 'BglII', site: 'AGATCT', cutTop: 1, cutBottom: 5 },
  { name: 'HindIII', site: 'AAGCTT', cutTop: 1, cutBottom: 5 },
  { name: 'XhoI', site: 'CTCGAG', cutTop: 1, cutBottom: 5 },
  { name: 'SalI', site: 'GTCGAC', cutTop: 1, cutBottom: 5 },
  { name: 'XbaI', site: 'TCTAGA', cutTop: 1, cutBottom: 5 },
  { name: 'SpeI', site: 'ACTAGT', cutTop: 1, cutBottom: 5 },
  { name: 'NheI', site: 'GCTAGC', cutTop: 1, cutBottom: 5 },
  { name: 'AflII', site: 'CTTAAG', cutTop: 1, cutBottom: 5 },
  { name: 'MluI', site: 'ACGCGT', cutTop: 1, cutBottom: 5 },
  { name: 'NcoI', site: 'CCATGG', cutTop: 1, cutBottom: 5 },
  { name: 'NotI', site: 'GCGGCCGC', cutTop: 2, cutBottom: 6 },
  { name: 'NdeI', site: 'CATATG', cutTop: 2, cutBottom: 4 },
  { name: 'ClaI', site: 'ATCGAT', cutTop: 2, cutBottom: 4 },
  // 3' overhangs.
  { name: 'KpnI', site: 'GGTACC', cutTop: 5, cutBottom: 1 },
  { name: 'SacI', site: 'GAGCTC', cutTop: 5, cutBottom: 1 },
  { name: 'PstI', site: 'CTGCAG', cutTop: 5, cutBottom: 1 },
  { name: 'SphI', site: 'GCATGC', cutTop: 5, cutBottom: 1 },
  // Blunt.
  { name: 'EcoRV', site: 'GATATC', cutTop: 3, cutBottom: 3 },
  { name: 'SmaI', site: 'CCCGGG', cutTop: 3, cutBottom: 3 },
  { name: 'PvuII', site: 'CAGCTG', cutTop: 3, cutBottom: 3 },
  { name: 'ScaI', site: 'AGTACT', cutTop: 3, cutBottom: 3 },
  { name: 'StuI', site: 'AGGCCT', cutTop: 3, cutBottom: 3 },
  { name: 'HpaI', site: 'GTTAAC', cutTop: 3, cutBottom: 3 },
  { name: 'DraI', site: 'TTTAAA', cutTop: 3, cutBottom: 3 },
];
