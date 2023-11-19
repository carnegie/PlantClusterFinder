BEGIN{
  # The script runs through a fasta file and writes a tag
  # "@PASCAL@" after every header line and every last
  # Sequence line of a contig/chromosome.
  # Pascal Schlapfer, January 19th 2016
  # Input: none (takes a file as input)
  # Output: writing a new file to standard output.
 
	currentmode = 0; # Identifies whether the last line
	                 # that was written, was a non-header
					 # line.
}
{
    
    firstChar = substr($0,1,1); # Reading in the first 
							    # Symbol of a line
    if ($firstChar ~ ">" ) { # If the line is a header line
		if (currentmode == 1) { # If the previous line was a
		                        # non-header line, then add
								# "@PASCAL@" to the sequence
								# of Nucleotides, and print it
								# to standard output.
			print vString, "@PASCAL@" # Add "@PASCAL@" to the
			                          # sequence of Nucleotides,
									  # and print it to standard
									  # output.
		}
		print $0, "@PASCAL@"; # Print the headerline with a
		                      # "@PASCAL@" tag.
		currentmode = 0; # Set the mode to be header.
    }
    else {
		if (currentmode == 0) { # If the previous line was a
		                        # headerline, then
			vString = $0; # Save the line for later use.
			currentmode = 1; # Set the mode to be non-header
		}
		else if (currentmode == 1) { # If the previous line was
		                             # not a header line, then
			print vString # Print the previous line that was stored,
			              # without the "@PASCAL@" tag.
			vString = $0; # Save the current line for later use.
		}
    }
}
END{
    if (currentmode == 1){ # If the very last line of the file was
	                       # a sequence line, then
		print vString, "@PASCAL@" # Print the line with a "@PASCAL@"
		                          # tag.
    }
}
