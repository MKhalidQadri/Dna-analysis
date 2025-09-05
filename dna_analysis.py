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

def longest_orf_in_sequence(records, seq_id, frame=1):
    """
    Find the longest ORF in a specific sequence by ID.
    """
    for record in records:
        if record.id == seq_id:
            orfs = find_orfs(record.seq, frame)
            if not orfs:
                return None
            longest_orf = max(orfs, key=lambda x: len(x[2]))
            return longest_orf
    return None

def find_repeats(records, n):
    """
    Find and count all repeats of length n in the sequences.
    """
    repeats = {}
    for record in records:
        seq = str(record.seq)
        for i in range(len(seq) - n + 1):
            repeat = seq[i:i+n]
            repeats[repeat] = repeats.get(repeat, 0) + 1
    max_repeat = max(repeats.items(), key=lambda x: x[1]) if repeats else (None, 0)
    return repeats, max_repeat

def print_separator():
    print("\n" + "="*60 + "\n")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python dna_analysis.py <multi-fasta-file>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    num_records, records = count_records(fasta_file)
    print_separator()
    print(f"Total Number of Records: {num_records}")
    print_separator()

    lengths, max_len, min_len, longest, shortest = sequence_lengths(records)
    print("Sequence Lengths:")
    for seq_id in sorted(lengths, key=lengths.get, reverse=True):
        print(f"  {seq_id}: {lengths[seq_id]} bp")
    print_separator()
    print(f"Longest Sequence Length: {max_len} bp")
    print(f"  Sequence(s): {', '.join(longest)}")
    print(f"Shortest Sequence Length: {min_len} bp")
    print(f"  Sequence(s): {', '.join(shortest)}")
    print_separator()


    frame = 1  # example frame
    longest_orf_len, longest_orf_id, longest_orf_seq, longest_orf_start = longest_orf_in_file(records, frame)
    print(f"Longest ORF in File (Frame {frame}):")
    print(f"  Length: {longest_orf_len} bp")
    print(f"  Sequence ID: {longest_orf_id}")
    print(f"  Start Position: {longest_orf_start}")
    print(f"  ORF Sequence: {longest_orf_seq}")
    print_separator()

    # Example: longest ORF in a given sequence ID
    if longest_orf_id:
        longest_orf = longest_orf_in_sequence(records, longest_orf_id, frame)
        if longest_orf:
            start, end, orf_seq = longest_orf
            print(f"Longest ORF in Sequence '{longest_orf_id}':")
            print(f"  Start: {start}, End: {end}")
            print(f"  ORF Sequence: {orf_seq}")
            print_separator()

    # Example: find repeats of length n
    n = 3  # example repeat length
    repeats, max_repeat = find_repeats(records, n)
    print(f"Repeats of Length {n}:")
    sorted_repeats = sorted(repeats.items(), key=lambda x: x[1], reverse=True)
    for repeat, count in sorted_repeats[:10]:  # Show top 10 repeats
        print(f"  {repeat}: {count} times")
    print_separator()
    print(f"Most Frequent Repeat of Length {n}: {max_repeat[0]} ({max_repeat[1]} times)")
    print_separator()




    # ...existing code...

frame2 = 2
longest_orf_len2, longest_orf_id2, longest_orf_seq2, longest_orf_start2 = longest_orf_in_file(records, frame2)
print(f"Longest ORF in File (Frame {frame2}):")
print(f"  Length: {longest_orf_len2} bp")
print(f"  Sequence ID: {longest_orf_id2}")
print(f"  Start Position: {longest_orf_start2}")
print(f"  ORF Sequence: {longest_orf_seq2}")
print_separator()

# ...existing code...

frame3 = 3
longest_orf_len3, longest_orf_id3, longest_orf_seq3, longest_orf_start3 = longest_orf_in_file(records, frame3)
print(f"Longest ORF in File (Frame {frame3}):")
print(f"  Length: {longest_orf_len3} bp")
print(f"  Sequence ID: {longest_orf_id3}")
print(f"  Start Position: {longest_orf_start3}")
print(f"  ORF Sequence: {longest_orf_seq3}")
print_separator()


# Find and print the most frequently occurring repeat of length 6 in all sequences
n = 6
repeats, max_repeat = find_repeats(records, n)
print(f"Most Frequent Repeat of Length {n}: {max_repeat[0]} ({max_repeat[1]} times)")
print_separator()


# Find and print all repeats of length 12 in all sequences
n = 12
repeats_12, max_repeat_12 = find_repeats(records, n)
print(f"Total different 12-base sequences: {len(repeats_12)}")
print(f"Most Frequent 12-base Repeat: {max_repeat_12[0]} (copies: {max_repeat_12[1]})")
print_separator()


# Find and print the most frequently occurring repeat of length 7 in all sequences
n = 7
repeats_7, max_repeat_7 = find_repeats(records, n)
print(f"Most Frequent Repeat of Length {n}: {max_repeat_7[0]} ({max_repeat_7[1]} times)")
print_separator()


# Find and print the longest forward ORF in the specified sequence
seq_id = "gi|142022655|gb|EQ086233.1|16"
longest_orf = None
longest_len = 0
longest_frame = None

for frame in [1, 2, 3]:
    orf = longest_orf_in_sequence(records, seq_id, frame)
    if orf:
        start, end, orf_seq = orf
        if len(orf_seq) > longest_len:
            longest_len = len(orf_seq)
            longest_orf = orf
            longest_frame = frame

if longest_orf:
    start, end, orf_seq = longest_orf
    print(f"Longest forward ORF in sequence '{seq_id}':")
    print(f"  Frame: {longest_frame}")
    print(f"  Start: {start}, End: {end}")
    print(f"  Length: {len(orf_seq)} bp")
    print(f"  ORF Sequence: {orf_seq}")
    print_separator()
else:
    print(f"No ORF found in sequence '{seq_id}'.")
    print_separator()



    # Find and print all repeats of length 12 in all sequences, and how many different 12-base sequences occur Max times
n = 12
repeats_12, max_repeat_12 = find_repeats(records, n)
Max = max_repeat_12[1]
most_frequent_12mers = [seq for seq, count in repeats_12.items() if count == Max]

print(f"Total different 12-base sequences: {len(repeats_12)}")
print(f"Most Frequent 12-base Repeat(s) (copies: {Max}):")
for seq in most_frequent_12mers:
    print(f"  {seq}")
print(f"Number of different 12-base sequences that occur {Max} times: {len(most_frequent_12mers)}")
print_separator()