/** Column-wise gap cleanup on a set of aligned sequences (a gap column is one
 * where every sequence has '-'). */

/**
 * Removes leading and trailing gap columns from an array of aligned sequences.
 * A column is trimmed only if ALL sequences have a gap (or end) there.
 * @param sequences Array of aligned sequences
 * @returns Array of trimmed sequences
 */
export function removeLeadingAndTrailingGapColumns(sequences: string[]): string[] {
  if (sequences.length === 0) return sequences;

  const maxLength = Math.max(...sequences.map(s => s.length));

  // Find first column with at least one non-gap character
  let firstColumn = maxLength;
  for (let col = 0; col < maxLength; col++) {
    const hasNonGap = sequences.some(seq => seq[col] && seq[col] !== '-');
    if (hasNonGap) {
      firstColumn = col;
      break;
    }
  }

  // Find last column with at least one non-gap character
  let lastColumn = -1;
  for (let col = maxLength - 1; col >= 0; col--) {
    const hasNonGap = sequences.some(seq => seq[col] && seq[col] !== '-');
    if (hasNonGap) {
      lastColumn = col;
      break;
    }
  }

  // If no non-gap characters found, return empty strings
  if (firstColumn >= maxLength || lastColumn < 0) {
    return sequences.map(() => '');
  }

  // Trim all sequences
  return sequences.map(seq => seq.substring(firstColumn, lastColumn + 1));
}

/**
 * Removes all gap columns from an array of aligned sequences.
 * A column is removed only if ALL sequences have a gap there.
 * @param sequences Array of aligned sequences
 * @returns Array of sequences with gap columns removed
 */
export function removeColumnsOfGaps(sequences: string[]): string[] {
  if (sequences.length === 0) return sequences;

  const maxLength = Math.max(...sequences.map(s => s.length));
  const columnsToKeep: number[] = [];

  // Identify columns to keep (at least one sequence has a non-gap character)
  for (let col = 0; col < maxLength; col++) {
    const hasNonGap = sequences.some(seq => seq[col] && seq[col] !== '-');
    if (hasNonGap) {
      columnsToKeep.push(col);
    }
  }

  // Build result sequences by keeping only the selected columns.
  return sequences.map(seq => columnsToKeep.map(col => seq[col] || '-').join(''));
}
