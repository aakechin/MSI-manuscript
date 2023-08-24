for read in bamfile.fetch(ref, start=start, end=end):
    if  not read.has_tag('MD'):
        continue
    readNum+=1
    readPos=read.pos
    # Stores position in read sequence
    x=0
    # Stores position in reference sequence that correspond to read
    y=0
    # We do not process unmapped or unaligned reads
    if read.cigartuples is None:
        continue
    # Now we need to find region with MSI-repeat
    # Because we got whole read
    # Go through all CIGAR regions
    for cigartuple in read.cigartuples:
        # If current region is match or deletion
        if cigartuple[0]==0 or cigartuple[0]==2:
            readPos+=cigartuple[1]
            y+=cigartuple[1]
        # If it is insertion, match or soft-clipped
        if cigartuple[0]==1 or cigartuple[0]==0 or cigartuple[0]==4:
            # We increase only position in read sequence
            x+=cigartuple[1]
        # If read position is more than start and
        # less than end
        if readPos>start and (readPos-1)<=end:
            # If there is deletion
            if cigartuple[0]==2:
                # We save it
                deletions.append(cigartuple[1])
        # If read position is equal to start
        # and first nucleotide has insertion
        if readPos==start and cigartuple[0]==1:
            # Extract reference sequence corresponding to read
            q=read.get_reference_sequence()
            # If previous nucleotide is the same as read's one
            # Else it can be the continue of insertion before MSI-region
            if q[y]==read.query_sequence[x]:
                # We save insertion
                insertions.append(cigartuple[1])
