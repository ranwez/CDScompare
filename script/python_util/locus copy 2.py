from array import array
class Locus:
    """
    Represents an annotation's locus identified from a GFF file.
    Uses __slots__ to reduce memory footprint.
    """
    __slots__ = ('name', 'mRNAs', 'start', 'end', 'phases')
    
    def __init__(self, name="", mRNAs=None, start=-1, end=-1, direction="", phases=None):
        """
        Initialize a Locus object.
        
        Args:
            name: ID of the locus
            mRNAs: Dictionary of mRNAs (ID -> list of CDS coordinates)
            start: Start coordinate of the locus
            end: End coordinate of the locus
            direction: Strand direction ('direct' or 'reverse')
            phases: Dictionary of phases (mRNA ID -> phase)
        """
        # Use the provided dictionary directly instead of copying
        self.name = name
        self.mRNAs = mRNAs if mRNAs is not None else {}
        self.start = start
        self.end = end
        #self.direction = direction

        # Only create a new dictionary if needed
        if phases is None:
            self.phases = {mRNA_id: 0 for mRNA_id in self.mRNAs}
        else:
            self.phases = phases  # Use the provided dictionary directly
    
    def get_mRNAs(self):
        """Return the mRNAs dictionary."""
        return self.mRNAs
    
    def set_mRNAs(self, value):
        """Set the mRNAs dictionary."""
        self.mRNAs = value
    
    def reverse(self, cluster_end):
        """
        Reverse the coordinates of all mRNAs relative to cluster_end.
        Used for loci on the reverse strand.
        
        Args:
            cluster_end: End position of the parent cluster
        """
        # Create a new dictionary to avoid modifying while iterating
        new_mRNAs = {}
        for mRNA_id, mRNA in self.mRNAs.items():
            # Using a list comprehension for better performance
            reversed_coords = [cluster_end - pos for pos in reversed(mRNA)]
            new_mRNAs[mRNA_id] = array('L', reversed_coords) 
        self.mRNAs = new_mRNAs
    
    def __str__(self):
        """Return a string representation of the locus."""
        return f"Locus '{self.name}' ({self.start}-{self.end})"
    
    def show_init(self):
        """Return a detailed string representation for debugging."""
        return f"Locus(name='{self.name}', mRNAs={self.mRNAs}, start={self.start}, end={self.end}')"