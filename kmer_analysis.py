"""
kmer_analysis.py

A command-line tool for counting k-mer frequencies and their subsequent
characters in DNA sequence files. This is a starting point for genome
assembly from sequenced fragments.

A k-mer is a substring of length k extracted from a DNA sequence. For each
k-mer, this tool records the total number of times it appears across all
sequences, as well as the frequency of each character that immediately
follows it.

Usage:
    python kmer_analysis.py <sequence_file> <k> <output_file>

Arguments:
    sequence_file : path to a text file with one DNA sequence per line
    k             : length of k-mers to count (integer)
    output_file   : path to write results to

Example:
    python kmer_analysis.py sequences.txt 2 results.txt
"""

import sys


def validate_sequence(sequence, k):
    """
    Check whether a DNA sequence is valid and long enough to produce k-mers.

    A valid sequence must:
    1. Be longer than k (needs at least one character after the last k-mer)
    2. Contain only valid DNA characters: A, T, G, C (upper or lowercase)
    3. Not contain digits, spaces, punctuation, or any other characters

    Parameters:

    sequence : str
        The DNA sequence string to validate.
    k : int
        The k-mer length. The sequence must be longer than k.

    Returns:

    bool
        True if the sequence is valid and long enough, False otherwise.

    Examples:

    >>> validate_sequence("ATGTCT", 2)
    True
    >>> validate_sequence("AT", 5)
    False
    >>> validate_sequence("ATG1CT", 2)
    False
    """
    # Reject sequences that are too short to produce any k-mer with a
    # following character. We need at least k+1 characters.
    if len(sequence) <= k:
        return False

    # Define the set of valid DNA nucleotide characters (upper and lowercase)
    valid_chars = set("ACGTacgt")

    # Check every character in the sequence against the valid set
    for nucleotide in sequence:
        if nucleotide not in valid_chars:
            return False

    return True


def update_kmer_count(kmer_data, kmer, next_char):
    """
    Add one occurrence of a k-mer and its following character to the data store.

    If the k-mer has not been seen before, it is initialized with a count of 1.
    If it has been seen before, its count is incremented by 1. The frequency of
    the character that follows the k-mer is also tracked separately.

    Parameters:

    kmer_data : dict
        A dictionary mapping each k-mer string to a sub-dictionary containing:
            - 'count'      : int, total number of times this k-mer has appeared
            - 'next_chars' : dict mapping each following character to its frequency
    kmer : str
        The k-mer string to record.
    next_char : str
        The single character that immediately follows this k-mer in the sequence.

    Returns:

    dict
        The updated kmer_data dictionary.

    Examples:

    >>> kmer_data = {}
    >>> kmer_data = update_kmer_count(kmer_data, "AT", "G")
    >>> kmer_data["AT"]["count"]
    1
    >>> kmer_data = update_kmer_count(kmer_data, "AT", "G")
    >>> kmer_data["AT"]["count"]
    2
    """
    # If this k-mer hasn't been seen before, initialize its entry.
    # Start count at 0 so the increment below brings it to 1 on first sight.
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}

    # Increment the total count for this k-mer
    kmer_data[kmer]['count'] += 1

    # If this following character hasn't been seen after this k-mer, initialize it
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0

    # Increment the count for this following character
    kmer_data[kmer]['next_chars'][next_char] += 1

    return kmer_data


def count_kmers_with_context(sequence, k):
    """
    Extract all k-mers and their following characters from a single sequence.

    Slides a window of length k across the sequence one position at a time.
    For each position i, the k-mer is sequence[i:i+k] and the following
    character is sequence[i+k]. The last k characters of the sequence are
    not counted as k-mers because they have no following character.

    Parameters:

    sequence : str
        A DNA sequence string. Should be validated before calling this function.
    k : int
        The length of k-mers to extract.

    Returns:

    dict
        A kmer_data dictionary mapping each k-mer to its count and next_chars.
        See update_kmer_count for the structure of this dictionary.

    Examples:

    >>> result = count_kmers_with_context("ATGT", 2)
    >>> result["AT"]["count"]
    1
    >>> result["AT"]["next_chars"]["G"]
    1
    """
    # Initialize an empty data store for this sequence
    kmer_data = {}

    # Slide across the sequence, stopping k positions before the end
    # so that every k-mer has a valid following character
    for i in range(len(sequence) - k):
        # Extract the k-mer starting at position i
        kmer = sequence[i:i+k]

        # The character immediately after the k-mer
        next_char = sequence[i+k]

        # Record this k-mer and its following character
        kmer_data = update_kmer_count(kmer_data, kmer, next_char)

    return kmer_data


def write_results_to_file(kmer_data, output_filename):
    """
    Write k-mer counts and next-character frequencies to an output file.

    Each line of the output file contains one k-mer, its total count, and
    the frequency of each character that follows it. K-mers are written in
    alphabetical order. Next characters within each line are also sorted
    alphabetically.

    Output format per line:
        <kmer> <total_count> <char1>:<freq1> <char2>:<freq2> ...

    Example output line:
        AT 3 G:2 C:1

    Parameters:

    kmer_data : dict
        A kmer_data dictionary as returned by count_kmers_with_context.
    output_filename : str
        Path to the output file to write. Will be created or overwritten.

    Returns:

    None

    Examples:

    >>> kmer_data = {"AT": {"count": 2, "next_chars": {"G": 2}}}
    >>> write_results_to_file(kmer_data, "output.txt")
    # output.txt will contain: AT 2 G:2
    """
    # Sort k-mers alphabetically for consistent, readable output
    sorted_kmers = sorted(kmer_data.keys())

    with open(output_filename, 'w') as f:
        for kmer in sorted_kmers:
            # Get the total count and next character frequencies for this k-mer
            total_count = kmer_data[kmer]['count']
            next_chars  = kmer_data[kmer]['next_chars']

            # Build the next-character string, sorted alphabetically by character
            next_char_str = " ".join(
                f"{char}:{freq}"
                for char, freq in sorted(next_chars.items())
            )

            # Write the k-mer, its total count, and next-char frequencies
            f.write(f"{kmer} {total_count} {next_char_str}\n")


def main():
    """
    Entry point for the k-mer analysis command-line tool.

    Reads DNA sequences from a file (one per line), counts k-mers and their
    following characters across ALL sequences, then writes the combined results
    to an output file.

    Command-line arguments (via sys.argv):
        1. sequence_file : path to input file with one sequence per line
        2. k             : k-mer length (integer)
        3. output_file   : path to write results

    Sequences that are invalid (contain non-DNA characters or are too short)
    are skipped with a warning message.

    Returns:

    None
    """
    # Read command-line arguments
    sequence_file = sys.argv[1]
    k             = int(sys.argv[2])
    output_file   = sys.argv[3]

    print(f"Reading sequences from {sequence_file}...")

    # Initialize a single shared kmer_data dictionary to accumulate counts
    # across ALL sequences in the file — not just the last one
    combined_kmer_data = {}

    with open(sequence_file, 'r') as f:
        for sequence in f:
            # Remove leading/trailing whitespace and newline characters
            sequence = sequence.strip()

            # Skip empty lines
            if not sequence:
                continue

            # Validate the sequence before processing
            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping invalid or too-short sequence: {sequence}")
                continue

            # Count k-mers for this sequence
            seq_kmer_data = count_kmers_with_context(sequence, k)

            # Merge this sequence's k-mer counts into the combined dictionary
            for kmer, data in seq_kmer_data.items():
                for next_char, freq in data['next_chars'].items():
                    # Add each occurrence individually to correctly update counts
                    for _ in range(freq):
                        combined_kmer_data = update_kmer_count(
                            combined_kmer_data, kmer, next_char
                        )

    # Write the combined results from all sequences to the output file
    write_results_to_file(combined_kmer_data, output_file)
    print(f"Results written to {output_file}")


if __name__ == '__main__':
    main()
