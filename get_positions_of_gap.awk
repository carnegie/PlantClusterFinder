{
  # Takes a FASTA file with header lines (>) followed by a single-line
  # sequence per chromosome/contig.
  # Input: FASTA from stdin or file
  # Output: contig, type (N or A), start, length

    # Non-header lines: process sequence
    if ( $0 !~ /^>/ ) {

                # Normalize bases:
                #   - All A/T/C/G (upper/lower) and '*' -> 'A'
                #   - All n/N -> 'N'
                seq = $0
                gsub(/[aAtTcCgG*]/, "A", seq)
                gsub(/[nN]/, "N", seq)

                len = length(seq)

                # Arrays to store N segments and A segments for this sequence
                n_cnt = 0
                a_cnt = 0

                # Single pass over the sequence to find contiguous N and A runs
                i = 1
                while (i <= len) {
                        base = substr(seq, i, 1)

                        if (base == "N") {
                                # Scan one contiguous N run
                                start = i
                                while (i <= len && substr(seq, i, 1) == "N") {
                                        i++
                                }
                                seglen = i - start

                                # Store this N run
                                n_cnt++
                                n_start[n_cnt] = start
                                n_len[n_cnt]   = seglen

                        } else if (base == "A") {
                                # Scan one contiguous A run (non-N bases)
                                start = i
                                while (i <= len && substr(seq, i, 1) == "A") {
                                        i++
                                }
                                seglen = i - start

                                # Store this A run
                                a_cnt++
                                a_start[a_cnt] = start
                                a_len[a_cnt]   = seglen

                        } else {
                                # Should not occur after normalization; skip defensively
                                i++
                        }
                }

                # First output all N segments (same order as original N-pass)
                for (k = 1; k <= n_cnt; k++) {
                        print vheader[1], "N", n_start[k], n_len[k]
                }

                # Then output all A segments (same order as original A-pass)
                for (k = 1; k <= a_cnt; k++) {
                        print vheader[1], "A", a_start[k], a_len[k]
                }

    } else {
                # Header line: store contig/chromosome ID in vheader[1]
                split($0, vheader, " ")
                gsub(/>/, "", vheader[1])
    }
}
