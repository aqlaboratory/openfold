from typing import Tuple,Iterable, Sequence,Optional

def parse_fasta(fasta_string: str) -> Tuple[Sequence[str], Sequence[str]]:
    """Parses FASTA string and returns list of strings with amino-acid sequences.
    Arguments:
        fasta_string: The string contents of a FASTA file.
    Returns:
        A tuple of two lists:
        * A list of sequences.
        * A list of sequence descriptions taken from the comment lines. In the
            same order as the sequences.
    """
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith(">"):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append("")
            continue
        elif line.startswith("#"):
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line

    return sequences, descriptions

def _convert_sto_seq_to_a3m(
    query_non_gaps: Sequence[bool], sto_seq: str
) -> Iterable[str]:
    for is_query_res_non_gap, sequence_res in zip(query_non_gaps, sto_seq):
        if is_query_res_non_gap:
            yield sequence_res
        elif sequence_res != "-":
            yield sequence_res.lower()

def convert_stockholm_to_a3m(
    stockholm_format: str, max_sequences: Optional[int] = None
) -> str:
    """Converts MSA in Stockholm format to the A3M format."""
    descriptions = {}
    sequences = {}
    reached_max_sequences = False

    for line in stockholm_format.splitlines():
        reached_max_sequences = (
            max_sequences and len(sequences) >= max_sequences
        )
        if line.strip() and not line.startswith(("#", "//")):
            # Ignore blank lines, markup and end symbols - remainder are alignment
            # sequence parts.
            seqname, aligned_seq = line.split(maxsplit=1)
            if seqname not in sequences:
                if reached_max_sequences:
                    continue
                sequences[seqname] = ""
            sequences[seqname] += aligned_seq

    for line in stockholm_format.splitlines():
        if line[:4] == "#=GS":
            # Description row - example format is:
            # #=GS UniRef90_Q9H5Z4/4-78            DE [subseq from] cDNA: FLJ22755 ...
            columns = line.split(maxsplit=3)
            seqname, feature = columns[1:3]
            value = columns[3] if len(columns) == 4 else ""
            if feature != "DE":
                continue
            if reached_max_sequences and seqname not in sequences:
                continue
            descriptions[seqname] = value
            if len(descriptions) == len(sequences):
                break

    # Convert sto format to a3m line by line
    a3m_sequences = {}
    # query_sequence is assumed to be the first sequence
    query_sequence = next(iter(sequences.values()))
    query_non_gaps = [res != "-" for res in query_sequence]
    for seqname, sto_sequence in sequences.items():
        a3m_sequences[seqname] = "".join(
            _convert_sto_seq_to_a3m(query_non_gaps, sto_sequence)
        )

    fasta_chunks = (
        f">{k} {descriptions.get(k, '')}\n{a3m_sequences[k]}"
        for k in a3m_sequences
    )
    return "\n".join(fasta_chunks) + "\n"  # Include terminating newline.