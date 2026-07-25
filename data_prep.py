def analyze_fasta(file_path):
    sequence_count = 0
    total_length = 0
    min_length = float('inf')
    max_length = 0

    sequence_ids = []
    current_seq_length = 0
    current_sequence = ""

    valid_fasta = True
    valid_nucleotides = True
    invalid_characters = set()

    # Valid nucleotide characters
    # Standard IUPAC nucleotide codes
    valid_chars = set("ATCGNRMWSYKVHDB")

    # Read FASTA file
    with open(file_path, 'r') as file:

        for line_number, line in enumerate(file, start=1):
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Check FASTA header
            if line.startswith('>'):

                # Save previous sequence statistics
                if current_seq_length > 0:
                    sequence_count += 1
                    total_length += current_seq_length

                    min_length = min(min_length, current_seq_length)
                    max_length = max(max_length, current_seq_length)

                # Get sequence identifier
                header = line[1:].strip()

                if not header:
                    valid_fasta = False
                    print(f"Warning: Empty FASTA header at line {line_number}")

                else:
                    # Get accession ID (first part of header)
                    accession_id = header.split()[0]
                    sequence_ids.append(accession_id)

                # Reset for new sequence
                current_seq_length = 0
                current_sequence = ""

            else:
                # Check if sequence appears before any FASTA header
                if sequence_count == 0 and not sequence_ids:
                    valid_fasta = False
                    print(
                        f"Warning: Sequence data found before FASTA header "
                        f"at line {line_number}"
                    )

                # Convert sequence to uppercase
                sequence = line.upper()

                # Check nucleotide characters
                for char in sequence:
                    if char not in valid_chars:
                        valid_nucleotides = False
                        invalid_characters.add(char)

                # Add sequence length
                current_seq_length += len(sequence)
                current_sequence += sequence

        # Save the final sequence
        if current_seq_length > 0:
            sequence_count += 1
            total_length += current_seq_length

            min_length = min(min_length, current_seq_length)
            max_length = max(max_length, current_seq_length)

    # Calculate average length
    avg_length = (
        total_length / sequence_count
        if sequence_count > 0
        else 0
    )

    # Check duplicate sequence IDs
    unique_ids = set(sequence_ids)
    duplicate_ids = []

    if len(sequence_ids) != len(unique_ids):
        seen = set()

        for seq_id in sequence_ids:
            if seq_id in seen and seq_id not in duplicate_ids:
                duplicate_ids.append(seq_id)
            seen.add(seq_id)

    # Print Task 2 Report
    print("\n--- TASK 2: DATA PREPARATION REPORT ---")

    print(f"\nTotal Number of Sequences: {sequence_count}")
    print(f"Total Base Pairs (bp): {total_length}")
    print(f"Shortest Sequence: {min_length} bp")
    print(f"Longest Sequence: {max_length} bp")
    print(f"Average Sequence Length: {avg_length:.1f} bp")

    # FASTA format validation
    print("\n--- FASTA FORMAT VALIDATION ---")

    if valid_fasta:
        print("FASTA Format: Valid")
    else:
        print("FASTA Format: Invalid")

    # Nucleotide validation
    print("\n--- NUCLEOTIDE VALIDATION ---")

    if valid_nucleotides:
        print("Nucleotide Characters: Valid (A, T, C, G, N only)")
    else:
        print("Nucleotide Characters: Invalid")
        print(f"Invalid Characters Found: {invalid_characters}")

    # Sequence ID validation
    print("\n--- SEQUENCE IDENTIFIER VALIDATION ---")

    print(f"Total Sequence IDs: {len(sequence_ids)}")
    print(f"Unique Sequence IDs: {len(unique_ids)}")

    if not duplicate_ids:
        print("Sequence Identifiers: All unique")
    else:
        print("Duplicate Sequence IDs Found:")
        for seq_id in duplicate_ids:
            print(f"- {seq_id}")

    # Complete genome length check
    print("\n--- COMPLETE GENOME LENGTH CHECK ---")

    # Approximate expected range for dengue virus complete genomes
    expected_min = 10000
    expected_max = 11000

    if min_length >= expected_min and max_length <= expected_max:
        print("Genome Length Check: All sequences are within the expected range.")
    else:
        print("Genome Length Check: Some sequences may require further inspection.")

    # Overall result
    print("\n--- TASK 2 SUMMARY ---")

    if (
        valid_fasta
        and valid_nucleotides
        and not duplicate_ids
        and min_length >= expected_min
        and max_length <= expected_max
    ):
        print("Data Preparation Status: PASSED")
    else:
        print("Data Preparation Status: REQUIRES REVIEW")


# Run the analysis
analyze_fasta("sequences (1).fasta")