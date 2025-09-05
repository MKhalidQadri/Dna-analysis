from Bio import SeqIO
from Bio.Seq import Seq

def count_records(fasta_file):
    """
    Return the number of records and a list of records in a FASTA file.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    return len(records), records

def sequence_lengths(records):
    """
    Return lengths of each sequence and info about longest and shortest.
    """
    lengths = {record.id: len(record.seq) for record in records}
    max_len = max(lengths.values())
    min_len = min(lengths.values())
    longest = [rid for rid, length in lengths.items() if length == max_len]
    shortest = [rid for rid, length in lengths.items() if length == min_len]
    return lengths, max_len, min_len, longest, shortest

def find_orfs(seq, frame=1):
    """
    Find all ORFs in a DNA sequence for a given frame.
    Returns a list of tuples: (start, end, orf_seq)
    """
    seq = seq[frame-1:]
    orfs = []
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    seq_len = len(seq)
    i = 0
    while i < seq_len - 2:
        codon = str(seq[i:i+3])
        if codon == start_codon:
            for j in range(i+3, seq_len-2, 3):
                stop_codon = str(seq[j:j+3])
                if stop_codon in stop_codons:
                    orfs.append((i+frame, j+3+frame-1, seq[i:j+3]))
                    i = j + 3
                    break
            else:
                i += 3
        else:
            i += 3
    return orfs

def longest_orf_in_file(records, frame=1):
    """
    Find the longest ORF in all records for a given frame.
    """
    longest_orf_len = 0
    longest_orf_seq = None
    longest_orf_id = None
    longest_orf_start = None
    for record in records:
        orfs = find_orfs(record.seq, frame)
        for start, end, orf_seq in orfs:
            orf_len = len(orf_seq)
            if orf_len > longest_orf_len:
                longest_orf_len = orf_len
                longest_orf_seq = orf_seq
                longest_orf_id = record.id
                longest_orf_start = start
    return longest_orf_len, longest_orf_id, longest_orf_seq, longest_orf_start

def find_repeats(records, n):
    """
    Find and count all repeats of length n in the sequences.
    Returns a dictionary of repeat:count and the most frequent repeat (repeat, count).
    """
    repeats = {}
    repeat_to_seqids = {}
    for record in records:
        seq = str(record.seq)
        for i in range(len(seq) - n + 1):
            repeat = seq[i:i+n]
            repeats[repeat] = repeats.get(repeat, 0) + 1
            if repeat not in repeat_to_seqids:
                repeat_to_seqids[repeat] = set()
            repeat_to_seqids[repeat].add(record.id)
    max_repeat = max(repeats.items(), key=lambda x: x[1]) if repeats else (None, 0)
    return repeats, max_repeat, repeat_to_seqids

def print_separator():
    print("\n" + "="*60 + "\n")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python Dnaanalysisupdated.py <multi-fasta-file>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    num_records, records = count_records(fasta_file)
    print_separator()
    print(f"1. Total Number of Records: {num_records}")
    print_separator()

    lengths, max_len, min_len, longest, shortest = sequence_lengths(records)
    print("2. Sequence Lengths:")
    for seq_id in sorted(lengths, key=lengths.get, reverse=True):
        print(f"  {seq_id}: {lengths[seq_id]} bp")
    print_separator()
    print(f"3. Longest Sequence Length: {max_len} bp")
    print(f"  Sequence(s): {', '.join(longest)}")
    print(f"4. Shortest Sequence Length: {min_len} bp")
    print(f"  Sequence(s): {', '.join(shortest)}")
    print_separator()

    # Longest ORF in file for each frame
    for idx, frame in enumerate([1, 2, 3], start=5):
        longest_orf_len, longest_orf_id, longest_orf_seq, longest_orf_start = longest_orf_in_file(records, frame)
        print(f"{idx}. Longest ORF in File (Frame {frame}):")
        print(f"  Length: {longest_orf_len} bp")
        print(f"  Sequence ID: {longest_orf_id}")
        print(f"  Start Position: {longest_orf_start}")
        print(f"  ORF Sequence: {longest_orf_seq}")
        print_separator()

    # 8. Longest ORF in File (All Frames)
    longest = (None, 0, None, None, None)  # (frame, length, id, start, seq)
    for frame in [1, 2, 3]:
        longest_orf_len, longest_orf_id, longest_orf_seq, longest_orf_start = longest_orf_in_file(records, frame)
        if longest_orf_len > longest[1]:
            longest = (frame, longest_orf_len, longest_orf_id, longest_orf_start, longest_orf_seq)
    frame, orf_len, orf_id, orf_start, orf_seq = longest
    print(f"8. Longest ORF in File (All Frames):")
    print(f"  Frame: {frame}")
    print(f"  Length: {orf_len} bp")
    print(f"  Sequence ID: {orf_id}")
    print(f"  Start Position: {orf_start}")
    print(f"  ORF Sequence: {orf_seq}")
    print_separator()

    # --- Combined Section: Repeats of different asked lengths (with numbering) ---
    repeat_lengths = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 18, 21, 24]
    for idx, n in enumerate(repeat_lengths, start=9):
        repeats_n, max_repeat_n, repeat_to_seqids = find_repeats(records, n)
        Max_n = max_repeat_n[1]
        most_frequent_nmers = [seq for seq, count in repeats_n.items() if count == Max_n]
        print(f"{idx}. Repeats of length {n}:")
        print(f"  Total different {n}-base sequences: {len(repeats_n)}")
        print(f"  Most Frequent Repeat(s) of Length {n} (copies: {Max_n}):")
        for seq in most_frequent_nmers:
            seq_names = ', '.join(sorted(repeat_to_seqids[seq]))
            print(f"    {seq} (found in: {seq_names})")
        print(f"  Number of different {n}-base sequences that occur {Max_n} times: {len(most_frequent_nmers)}")
        # Show how many sequences occur different times (frequency distribution)
        freq_dist = {}
        for count in repeats_n.values():
            freq_dist[count] = freq_dist.get(count, 0) + 1
        print("  Number of different base sequences that occur different times:")
        for count in sorted(freq_dist):
            print(f"    {freq_dist[count]} sequence(s) occur {count} times")
        print_separator()

