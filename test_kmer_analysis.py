# test_kmer_analysis.py
# Pytest test suite for kmer_analysis.py
#
# Tests are organized by function. Each function in kmer_analysis.py has
# multiple tests covering normal cases, edge cases, and expected failures.
#
# Run with:
#   pytest test_kmer_analysis.py
#   pytest test_kmer_analysis.py -v   (verbose output)

import pytest
import os
from kmer_analysis import (
    validate_sequence,
    update_kmer_count,
    count_kmers_with_context,
    write_results_to_file
)


# Tests for validate_sequence(sequence, k)

class TestValidateSequence:
    """Tests for the validate_sequence function."""

    def test_valid_sequence_returns_true(self):
        """A normal DNA sequence longer than k should return True."""
        assert validate_sequence("ATGTCTGTCTGAA", 2) == True

    def test_sequence_shorter_than_k_returns_false(self):
        """A sequence shorter than k cannot produce any k-mers."""
        assert validate_sequence("AT", 5) == False

    def test_sequence_equal_to_k_returns_false(self):
        """A sequence exactly equal to k has no character after the k-mer."""
        assert validate_sequence("AT", 2) == False

    def test_sequence_with_digits_returns_false(self):
        """Sequences containing digits are not valid DNA."""
        assert validate_sequence("ATG1CT", 2) == False

    def test_sequence_with_invalid_characters_returns_false(self):
        """Sequences containing non-DNA characters should be rejected."""
        assert validate_sequence("ATGXCT", 2) == False

    def test_sequence_with_spaces_returns_false(self):
        """Sequences containing spaces should be rejected."""
        assert validate_sequence("ATG CT", 2) == False

    def test_lowercase_sequence(self):
        """Lowercase DNA sequences should be handled."""
        assert validate_sequence("atgtct", 2) == True

    def test_empty_sequence_returns_false(self):
        """An empty sequence should return False."""
        assert validate_sequence("", 2) == False

    def test_valid_sequence_with_k_equals_1(self):
        """A sequence of length 2 with k=1 should be valid."""
        assert validate_sequence("AT", 1) == True

    def test_sequence_with_punctuation_returns_false(self):
        """Sequences containing punctuation should be rejected."""
        assert validate_sequence("ATG-CT", 2) == False


# Tests for update_kmer_count(kmer_data, kmer, next_char)

class TestUpdateKmerCount:
    """Tests for the update_kmer_count function."""

    def test_new_kmer_count_is_one(self):
        """A new kmer should have a count of exactly 1."""
        kmer_data = {}
        kmer_data = update_kmer_count(kmer_data, "AT", "G")
        assert kmer_data["AT"]["count"] == 1

    def test_existing_kmer_count_increments(self):
        """Seeing a kmer a second time should increment count to 2."""
        kmer_data = {}
        kmer_data = update_kmer_count(kmer_data, "AT", "G")
        kmer_data = update_kmer_count(kmer_data, "AT", "G")
        assert kmer_data["AT"]["count"] == 2

    def test_new_next_char_count_is_one(self):
        """A new next_char for a kmer should have frequency 1."""
        kmer_data = {}
        kmer_data = update_kmer_count(kmer_data, "AT", "G")
        assert kmer_data["AT"]["next_chars"]["G"] == 1

    def test_existing_next_char_increments(self):
        """Seeing the same next_char twice should increment its count to 2."""
        kmer_data = {}
        kmer_data = update_kmer_count(kmer_data, "AT", "G")
        kmer_data = update_kmer_count(kmer_data, "AT", "G")
        assert kmer_data["AT"]["next_chars"]["G"] == 2

    def test_multiple_next_chars_tracked_separately(self):
        """Different next_chars for the same kmer should be tracked independently."""
        kmer_data = {}
        kmer_data = update_kmer_count(kmer_data, "AT", "G")
        kmer_data = update_kmer_count(kmer_data, "AT", "C")
        assert kmer_data["AT"]["next_chars"]["G"] == 1
        assert kmer_data["AT"]["next_chars"]["C"] == 1

    def test_multiple_kmers_tracked_independently(self):
        """Different kmers should be tracked independently."""
        kmer_data = {}
        kmer_data = update_kmer_count(kmer_data, "AT", "G")
        kmer_data = update_kmer_count(kmer_data, "TG", "C")
        assert "AT" in kmer_data
        assert "TG" in kmer_data
        assert kmer_data["AT"]["count"] == 1
        assert kmer_data["TG"]["count"] == 1

    def test_returns_updated_dict(self):
        """Function should return the updated kmer_data dictionary."""
        kmer_data = {}
        result = update_kmer_count(kmer_data, "AT", "G")
        assert isinstance(result, dict)
        assert "AT" in result


# Tests for count_kmers_with_context(sequence, k)

class TestCountKmersWithContext:
    """Tests for the count_kmers_with_context function."""

    def test_basic_kmer_counting(self):
        """Simple sequence should produce correct kmer counts."""
        # ATGT with k=2: AT->G, TG->T
        result = count_kmers_with_context("ATGT", 2)
        assert "AT" in result
        assert "TG" in result

    def test_repeated_kmer_counted_correctly(self):
        """A kmer appearing twice should have count of 2."""
        # ATAT with k=2: AT->A, TA->T, AT (end - no next char so only 2 kmers)
        result = count_kmers_with_context("ATATAT", 2)
        assert result["AT"]["count"] == 2

    def test_next_char_recorded_correctly(self):
        """The character after each kmer should be recorded correctly."""
        result = count_kmers_with_context("ATGT", 2)
        assert "G" in result["AT"]["next_chars"]
        assert result["AT"]["next_chars"]["G"] == 1

    def test_k_equals_1(self):
        """Should work correctly with k=1."""
        result = count_kmers_with_context("ATGT", 1)
        assert "A" in result
        assert "T" in result
        assert "G" in result

    def test_last_kmer_excluded(self):
        """The last k characters have no following character and should not be counted."""
        # "ATG" with k=2: only "AT" has a next char (G), "TG" is the last kmer
        result = count_kmers_with_context("ATG", 2)
        assert "AT" in result
        assert "TG" not in result

    def test_known_example_from_assignment(self):
        """Test the example given in the assignment: ATGTCTGTCTGAA with k=2."""
        result = count_kmers_with_context("ATGTCTGTCTGAA", 2)
        # TG appears three times and is followed by T both times
        assert result["TG"]["count"] == 3
        assert result["TG"]["next_chars"]["T"] == 2
        assert result["TG"]["next_chars"]["A"] == 1
        
    def test_returns_dict(self):
        """Function should return a dictionary."""
        result = count_kmers_with_context("ATGT", 2)
        assert isinstance(result, dict)


# Tests for write_results_to_file(kmer_data, output_filename)

class TestWriteResultsToFile:
    """Tests for the write_results_to_file function."""

    def test_output_file_is_created(self, tmp_path):
        """The output file should be created after writing."""
        kmer_data = {"AT": {"count": 1, "next_chars": {"G": 1}}}
        output_file = str(tmp_path / "output.txt")
        write_results_to_file(kmer_data, output_file)
        assert os.path.exists(output_file)

    def test_kmer_appears_in_output(self, tmp_path):
        """Each kmer should appear in the output file."""
        kmer_data = {"AT": {"count": 2, "next_chars": {"G": 2}}}
        output_file = str(tmp_path / "output.txt")
        write_results_to_file(kmer_data, output_file)
        with open(output_file) as f:
            content = f.read()
        assert "AT" in content

    def test_total_count_in_output(self, tmp_path):
        """The total kmer count should appear in the output."""
        kmer_data = {"AT": {"count": 2, "next_chars": {"G": 2}}}
        output_file = str(tmp_path / "output.txt")
        write_results_to_file(kmer_data, output_file)
        with open(output_file) as f:
            content = f.read()
        assert "2" in content

    def test_next_char_frequency_in_output(self, tmp_path):
        """Next character frequencies should appear in the output."""
        kmer_data = {"AT": {"count": 1, "next_chars": {"G": 1}}}
        output_file = str(tmp_path / "output.txt")
        write_results_to_file(kmer_data, output_file)
        with open(output_file) as f:
            content = f.read()
        assert "G:1" in content

    def test_output_is_sorted_alphabetically(self, tmp_path):
        """Kmers should be written in alphabetical order."""
        kmer_data = {
            "TG": {"count": 1, "next_chars": {"T": 1}},
            "AT": {"count": 1, "next_chars": {"G": 1}}
        }
        output_file = str(tmp_path / "output.txt")
        write_results_to_file(kmer_data, output_file)
        with open(output_file) as f:
            lines = f.readlines()
        assert lines[0].startswith("AT")
        assert lines[1].startswith("TG")

    def test_multiple_next_chars_in_output(self, tmp_path):
        """Multiple next characters for a kmer should all appear in output."""
        kmer_data = {"AT": {"count": 2, "next_chars": {"G": 1, "C": 1}}}
        output_file = str(tmp_path / "output.txt")
        write_results_to_file(kmer_data, output_file)
        with open(output_file) as f:
            content = f.read()
        assert "G:1" in content
        assert "C:1" in content

    def test_empty_kmer_data_creates_empty_file(self, tmp_path):
        """Empty kmer_data should produce an empty output file."""
        output_file = str(tmp_path / "output.txt")
        write_results_to_file({}, output_file)
        with open(output_file) as f:
            content = f.read()
        assert content == ""


# Integration test

class TestIntegration:
    """End-to-end tests running the full pipeline."""

    def test_full_pipeline_single_sequence(self, tmp_path):
        """Full pipeline should produce correct output for a single sequence."""
        # Create input file
        input_file = str(tmp_path / "sequences.txt")
        with open(input_file, "w") as f:
            f.write("ATGTCTGTCTGAA\n")

        # Run counting
        result = count_kmers_with_context("ATGTCTGTCTGAA", 2)

        # Write output
        output_file = str(tmp_path / "output.txt")
        write_results_to_file(result, output_file)

        # Check output
        with open(output_file) as f:
            content = f.read()

        assert "TG" in content
        assert "AT" in content

    def test_full_pipeline_multiple_sequences(self, tmp_path):
        """Kmers should be accumulated correctly across multiple sequences."""
        from kmer_analysis import count_kmers_with_context, update_kmer_count

        sequences = ["ATGT", "ATGC"]
        combined = {}
        for seq in sequences:
            seq_data = count_kmers_with_context(seq, 2)
            for kmer, data in seq_data.items():
                for next_char, freq in data["next_chars"].items():
                    for _ in range(freq):
                        combined = update_kmer_count(combined, kmer, next_char)

        # AT appears in both sequences
        assert combined["AT"]["count"] == 2
