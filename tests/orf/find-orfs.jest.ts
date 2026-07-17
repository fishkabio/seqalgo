import { reverseComplement } from '../../src/sequence/sequence.utils';
import { findOrfs } from '../../src/orf/find-orfs';

describe('findOrfs', () => {
  test('finds a simple forward ORF from ATG to the stop, in genomic half-open coords', () => {
    //            0         1         2
    //            0123456789012345678901
    const seq = 'AAATGAAACCCGGGTAACC'; // ATG at index 2, TAA stop at index 14..17
    const orfs = findOrfs(seq, { minAminoAcids: 1 });
    const fwd = orfs.filter(o => o.strand === 1);
    expect(fwd).toHaveLength(1);
    expect(fwd[0]).toMatchObject({ start: 2, end: 17, strand: 1, hasStop: true, startCodon: 'ATG', stopCodon: 'TAA' });
    // ATG AAA CCC GGG = 4 codons before the stop.
    expect(fwd[0].aaLength).toBe(4);
    // The slice is a clean ATG…TAA.
    expect(seq.slice(fwd[0].start, fwd[0].end)).toBe('ATGAAACCCGGGTAA');
  });

  test('respects the minimum amino-acid length', () => {
    const seq = 'ATGAAATAA'; // ORF of 2 aa (ATG, AAA) then TAA
    expect(findOrfs(seq, { minAminoAcids: 2, includeReverse: false })).toHaveLength(1);
    expect(findOrfs(seq, { minAminoAcids: 3, includeReverse: false })).toHaveLength(0);
  });

  test('takes the FIRST start codon in a stop-delimited region (longest ORF for that stop)', () => {
    // ATG ... ATG ... TAA — the ORF should begin at the first ATG (index 0), not the inner one.
    const seq = 'ATGCCCATGGGGTAA';
    const [orf] = findOrfs(seq, { minAminoAcids: 1, includeReverse: false });
    expect(orf.start).toBe(0);
    expect(orf.aaLength).toBe(4); // ATG CCC ATG GGG
  });

  test('maps a reverse-strand ORF back to plus-strand coordinates', () => {
    // Build a sequence whose reverse complement contains a clean ORF, then check the mapping.
    const rcOrf = 'ATGTTTGGGTAA'; // a 3-aa ORF on the minus strand's 5'->3'
    const seq = reverseComplement(rcOrf); // put it on the minus strand of `seq`
    const minus = findOrfs(seq, { minAminoAcids: 1 }).filter(o => o.strand === -1);
    expect(minus).toHaveLength(1);
    expect(minus[0]).toMatchObject({ start: 0, end: seq.length, strand: -1, startCodon: 'ATG', stopCodon: 'TAA' });
    // Reading the minus strand 5'->3' over the reported span reproduces the ORF.
    expect(reverseComplement(seq.slice(minus[0].start, minus[0].end))).toBe(rcOrf);
  });

  test('finds ORFs in non-zero frames', () => {
    const seq = 'GG' + 'ATGAAATAA'; // ORF starts at index 2 (frame 2)
    const [orf] = findOrfs(seq, { minAminoAcids: 1, includeReverse: false });
    expect(orf).toMatchObject({ start: 2, frame: 2, strand: 1 });
  });

  test('honours a custom start-codon set and can skip the reverse strand', () => {
    const seq = 'GTGAAATAA'; // GTG start
    expect(findOrfs(seq, { minAminoAcids: 1, includeReverse: false })).toHaveLength(0); // ATG-only default
    const withGtg = findOrfs(seq, { minAminoAcids: 1, includeReverse: false, startCodons: ['ATG', 'GTG'] });
    expect(withGtg).toHaveLength(1);
    expect(withGtg[0].startCodon).toBe('GTG');
  });

  test('with requireStop:false, an unterminated ORF runs to the end of its frame', () => {
    const seq = 'ATGAAACCC'; // no stop
    expect(findOrfs(seq, { minAminoAcids: 1, includeReverse: false })).toHaveLength(0); // requireStop default
    const open = findOrfs(seq, { minAminoAcids: 1, includeReverse: false, requireStop: false });
    expect(open).toHaveLength(1);
    expect(open[0]).toMatchObject({ start: 0, end: 9, hasStop: false, stopCodon: '', aaLength: 3 });
  });

  test('reads RNA input as DNA (U→T)', () => {
    const orfs = findOrfs('AUGAAACCCUAA', { minAminoAcids: 1, includeReverse: false });
    expect(orfs).toHaveLength(1);
    expect(orfs[0]).toMatchObject({ start: 0, end: 12, startCodon: 'ATG', stopCodon: 'TAA', aaLength: 3 });
  });

  test('ignores whitespace in the input sequence', () => {
    const orfs = findOrfs('ATG AAA\nCCC GGG TAA', { minAminoAcids: 1, includeReverse: false });
    expect(orfs).toHaveLength(1);
    expect(orfs[0]).toMatchObject({ start: 0, end: 15, aaLength: 4 });
  });
});
