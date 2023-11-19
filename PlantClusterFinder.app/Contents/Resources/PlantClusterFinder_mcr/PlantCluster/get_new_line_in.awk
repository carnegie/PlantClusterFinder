{
  # Replaces all "@PASCAL@" tags with a new line character
  # and sends the line to the standard output.
  # Pascal Schlapfer, January 19th 2016
    gsub("@PASCAL@","\n") # Replace the tag.
    print $0 # Print to standard output.
}
