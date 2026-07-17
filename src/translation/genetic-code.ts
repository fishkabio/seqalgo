/** A genetic code: uppercase DNA codon -> single-letter amino acid, with `*` for a stop codon. */
export type GeneticCode = Readonly<Record<string, string>>;

/** One NCBI translation table: its `transl_table` id, official NCBI name, and codon->amino-acid map. */
export interface GeneticCodeTable {
  /** NCBI `transl_table` id (1 = Standard). Ids are not contiguous (7, 8, 17-20 are retired). */
  readonly id: number;
  /** Official NCBI name (whitespace-collapsed; may list several organism groups). */
  readonly name: string;
  readonly code: GeneticCode;
}

/**
 * The 64 codons in NCBI `ncbieaa` order — base 1 outermost, base 3 innermost, each cycling T,C,A,G.
 * This is the order every table's amino-acid string below is written in, so zipping the two yields
 * the codon->amino-acid map.
 */
const CODONS: readonly string[] = ((): string[] => {
  const bases = 'TCAG';
  const out: string[] = [];
  for (const b1 of bases) for (const b2 of bases) for (const b3 of bases) out.push(b1 + b2 + b3);
  return out;
})();

/** Build a {@link GeneticCode} map from an NCBI 64-char amino-acid string (one letter per codon). */
function buildCode(ncbieaa: string): GeneticCode {
  const code: Record<string, string> = {};
  for (let i = 0; i < CODONS.length; i++) code[CODONS[i]] = ncbieaa[i];
  return code;
}

/**
 * Raw NCBI genetic-code data — [id, name, ncbieaa] — taken verbatim from NCBI's `gc.prt`
 * (https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt). The `ncbieaa` string is the amino acid
 * for each of the 64 codons in {@link CODONS} order; `*` marks a stop. Start-codon data (`sncbieaa`)
 * is intentionally omitted — plain frame translation uses only the amino-acid table.
 */
const RAW: readonly (readonly [number, string, string])[] = [
  [1, 'Standard', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [2, 'Vertebrate Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG'],
  [3, 'Yeast Mitochondrial', 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [
    4,
    'Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma',
    'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
  ],
  [5, 'Invertebrate Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG'],
  [
    6,
    'Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear',
    'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
  ],
  [
    9,
    'Echinoderm Mitochondrial; Flatworm Mitochondrial',
    'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
  ],
  [10, 'Euplotid Nuclear', 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [11, 'Bacterial, Archaeal and Plant Plastid', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [12, 'Alternative Yeast Nuclear', 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [13, 'Ascidian Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG'],
  [14, 'Alternative Flatworm Mitochondrial', 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG'],
  [15, 'Blepharisma Macronuclear', 'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [16, 'Chlorophycean Mitochondrial', 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [21, 'Trematode Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG'],
  [22, 'Scenedesmus obliquus Mitochondrial', 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [23, 'Thraustochytrium Mitochondrial', 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [24, 'Rhabdopleuridae Mitochondrial', 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG'],
  [
    25,
    'Candidate Division SR1 and Gracilibacteria',
    'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
  ],
  [26, 'Pachysolen tannophilus Nuclear', 'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [27, 'Karyorelict Nuclear', 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [28, 'Condylostoma Nuclear', 'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [29, 'Mesodinium Nuclear', 'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [30, 'Peritrich Nuclear', 'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [31, 'Blastocrithidia Nuclear', 'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [32, 'Balanophoraceae Plastid', 'FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'],
  [33, 'Cephalodiscidae Mitochondrial', 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG'],
];

/** Every NCBI translation table, in id order — the source of truth for {@link getGeneticCode}. */
export const GENETIC_CODES: readonly GeneticCodeTable[] = RAW.map(([id, name, ncbieaa]) => ({
  id,
  name,
  code: buildCode(ncbieaa),
}));

const BY_ID: ReadonlyMap<number, GeneticCodeTable> = new Map(GENETIC_CODES.map(t => [t.id, t]));

/** The NCBI table with this `transl_table` id, or `undefined` if no such table exists. */
export function getGeneticCodeTable(id: number): GeneticCodeTable | undefined {
  return BY_ID.get(id);
}

/** The codon->amino-acid map for this `transl_table` id, or `undefined` if no such table exists. */
export function getGeneticCode(id: number): GeneticCode | undefined {
  return BY_ID.get(id)?.code;
}

/**
 * NCBI translation table 1 (the Standard Code). Also the codon->amino-acid map for table 11
 * (Bacterial/Plant Plastid), which only differs from the standard code in alternative start codons.
 */
export const STANDARD_GENETIC_CODE: GeneticCode = buildCode(
  'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
);
