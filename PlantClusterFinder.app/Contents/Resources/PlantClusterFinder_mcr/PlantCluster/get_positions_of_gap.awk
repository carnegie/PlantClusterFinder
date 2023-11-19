{
  # Takes a fasta-file as input that is organized by header
  # lines, followed by a single line of the whole sequence
  # of the chromosome/contig.
  # Pascal Schlapfer, January 19th 2016
  # Input: none (takes a file as input)
  # Output: writing a new file to standard output.
  
    if ( substr($0,1,1) !~ ">" ) { # If the line is not a
	                               # header line, then
		# Make sequences consisting of only "A", "T", "C",
		# "G", "a", "t", "c", "g", and "*" identifiable:
		gsub("a","A",$0); # replace all small letter "a" with
		                  # capital "A"s.
		gsub(/\*/,"A",$0); # the same with asterisk-character.
		gsub("c","A",$0); # the same with "c".
		gsub("g","A",$0); # the same with "g".
		gsub("t","A",$0); # the same with "t".
		gsub("C","A",$0); # the same with "C".
		gsub("G","A",$0); # the same with "G".
		gsub("T","A",$0); # the same with "T".
		# Make sequences consisting of only "N" and "n" 
		# identifiable.
		gsub("n","N",$0);# replace all small letter "n" with
		                 # capital "N"s.
		vS1 = $0; # Store the line to be able to use it 
				  # twice.
		# First, identify all N-sequence-gaps:
		vPlace = match(vS1,/N+/,N_arr); # Identify the first
		                                # occurence of an "N" in
										# the sequence and identify
										# the length of the
										# "N"-sequence".
		vRepl = 0; # Initialize the length of the already searched
		           # part of the line by setting it to zero.
		while (N_arr[0] != "") { # while during a search, an "N"-sequence-gap was
		                         # found, then
			# Print out the current contig we are investigating, an "N" that
			# identifies that this is a N-sequence-gap, the start of the gap,
			# and the length of the "N"-sequence-gap.
			print vheader[1], "N", (N_arr[0, "start"] + vRepl), N_arr[0, "length"];
			# From the line, cut off everything from the 1st character of the
			# line until the last character of the found "N"-sequence-gap.
			vS1 = substr(vS1,N_arr[0, "start"]+N_arr[0, "length"])
			# increase the length of the already searched part of the line
			# by the strech that was just cut away.
			vRepl = vRepl + N_arr[0, "start"] + N_arr[0, "length"]-1;
			# Find a new "N"-sequence word.
			vPlace = match(vS1,/N+/,N_arr); # vPlace is zero, if no "N"-sequence
											# gap was found.
		}
		# Second, identify all ATCGatcg*-sequences. This is the same as above,
		# with the only difference that not "N"-sequences are searched but
		# "A"-sequences.
		vS1 = $0; # Restore the line.
		vPlace = match(vS1,/A+/,A_arr);
		vRepl = 0;
		while (A_arr[0] != "") {
			print vheader[1], "A", (A_arr[0, "start"] + vRepl), A_arr[0, "length"];
			vS1 = substr(vS1,A_arr[0, "start"]+A_arr[0, "length"]);
			vRepl = vRepl + A_arr[0, "start"] + A_arr[0, "length"]-1;
			vPlace = match(vS1,/A+/,A_arr);
		} 
    }
    else { # If the line is a header line, then
		split($0,vheader," "); # Split off the first word, containing the description
							   # of the contig/chromosome and save it in vheader
		gsub(/>/,"",vheader[1]); # Replace the leading ">" character of the variant.
    }
}
	
