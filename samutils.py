import re
import sys


alphabet = ['A','C','G','T','N','R','Y']
alphabetString = "".join(alphabet)
complementMapping = {'A':'T','C':'G','T':'A','G':'C',
                     '-':'-','N':'N','Y':'Y','R':'R'}
#
# Utility functions to parse MD String
#  of an alignment
def getMDStringComps(mdString):
    regExMDInit = '([0-9]+)'
    regExMDRest =  '([' + alphabetString + ']+|\^[' + alphabetString + ']+)([0-9]+)'
    # Parse MD String to get its components
    reM = re.match(regExMDInit, mdString)
    mdComps = [reM.group(1)]
    comps = re.findall(regExMDRest, mdString)
    for c,d in comps:
        mdComps += [c,d]
    return mdComps


#
# Utility function to parse CIGAR string of an
#  alignment
def getCIGARPairs(cigarString):
    regExCIGAR = '([0-9]+)([MIDNSHPX=])'
    cigarPairs = []
    # Parse the CIGAR String to get the components
    cpairs = re.findall(regExCIGAR, cigarString)
    for l,m in cpairs:
        cigarPairs += [[int(l),m]]
    return cigarPairs

#
# Builds the reference alignment from
#  (i)  CIGAR Pairs (parsed from getCIGARPairs),
#  (ii) MD components (parsed from getMDStringComps(...))
# and the read string
def getRefAlignment(cigarPairs,mdComps,readString):
    # Fix up the ref alignment string
    sofar = 0
    refAlign = ''
    [fl,fc] = cigarPairs[0]
    if fc == 'S':
        sofar += fl
        cigarPairs = cigarPairs[1:]
    for [l,c] in cigarPairs:
        if c == 'M':
            refAlign += readString[sofar:sofar+l]
            sofar += l
        elif c == 'S':
            sofar += l
        elif c == 'I':
            refAlign += l * '-' # nothin in ref
            sofar += l
        elif c == 'D':
            refAlign += l * 'N' # something in ref
    # Update the ref align based on MD string
    sofar = 0
    alnidx = 0
    alnlen = len(refAlign)
    refAlign2 = ''
    #print refAlign
    for c in mdComps:
        if c.isdigit():
            # If we have digits, we just move along
            # copying the strings
            limit = sofar + int(c)
            while sofar < limit:
                refAlign2 += refAlign[alnidx]
                if refAlign[alnidx] != '-':
                    sofar += 1
                alnidx += 1
        elif c[0] == '^': #deleted in read
            # If deleted in read, it is present
            refAlign2 += c[1:]
            alnidx += len(c[1:])
        else: # Substitution
            refAlign2 += c
            alnidx += len(c)
        while alnidx < alnlen and refAlign[alnidx] == '-':
            refAlign2 += refAlign[alnidx]
            alnidx += 1
    #print refAlign2
    return refAlign2

#
# Get the reference alignment from
#   - genome/chromosome string
#   - CIGAR info (parsed from CIGAR string with getCIGARPairs(...)
#   - starting position at which the genome is aligned at
#   - read string
def getGenomeRefAlign(gstr,gpos,cigarPairs,readString):
    sofar = gpos-1
    refAlign = ''
    [fl,fc] = cigarPairs[0]
    if fc == 'S':
        cigarPairs = cigarPairs[1:]
    for [l,c] in cigarPairs:
        if c == 'M':
            refAlign += gstr[sofar:sofar+l]
            sofar += l
        elif c == 'S':
            sofar += l
            break # Skipping should be the last
        elif c == 'I':
            refAlign += l * '-' # nothin in ref
        elif c == 'D':
            refAlign += gstr[sofar:sofar+l] # something in ref
            sofar += l
    #x = len(refAlign)
    #print gstr[gpos-1:gpos+x]
    return refAlign

#
# Get the alignment of read string from
#  CIGAR info (parsed from CIGAR string with getCIGARPairs(..) fn),
#  and read string
def getReadAlignment(cigarPairs,readString):
    # Fix up the read alignment string
    sofar = 0
    readAlign = ''
    skiplength = 0
    readlen = len(readString)
    skipbegin = 0
    skipend = readlen
    [fl,fc] = cigarPairs[0]
    if fc == 'S':
        skipbegin = fl
        skiplength += fl
        sofar += fl
        cigarPairs = cigarPairs[1:]
    for [l,c] in cigarPairs:
        if c == 'M':
            readAlign += readString[sofar:sofar+l]
            sofar += l
        elif c == 'S':
            skipend = sofar
            sofar += l
            skiplength += l
        elif c == 'I':
            readAlign += readString[sofar:sofar+l] # some thing in read
            sofar += l
        elif c == 'D':
            readAlign += l * '-' # nothing in read
    #print readAlign
    #if skiplength > 0: #(0.1 * readlen):
    #    raise SkipError(str(skiplength))
    return (readAlign,skipbegin,skipend)

def getSAMAlignment(readString, cigarString, mdString,
                    gstr = None, gpos = None):
    if mdString is None and gstr is None:
        return [ None, None, None, None ]
    if gstr is not None and gpos is None:
        return [ None, None, None, None ]
    readLength = len(readString)
    # Parse CIGAR and MD String
    cigarPairs = getCIGARPairs(cigarString)
    # Fix up the read and ref alignment string
    (readAlign,skipbegin,skipend) = getReadAlignment(cigarPairs,readString)
    refAlign = ''
    if mdString != None:
        mdComps = getMDStringComps(mdString)
        refAlign = getRefAlignment(cigarPairs,mdComps,readString)
        #print mdString,mdComps,refAlign
    else:
        refAlign = getGenomeRefAlign(gstr,gpos,cigarPairs,readString)
    if(len(refAlign) != len(readAlign)):
        print 'length doesnt match'
        print readAlign
        print refAlign
    assert len(refAlign) == len(readAlign)
    #print readAlign
    #print refAlign
    return [refAlign, readAlign,skipbegin,skipend]
