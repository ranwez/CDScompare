# -*- coding: utf-8 -*-

from attrs import define, field

@define(slots=True, eq=True)
class MatchScore:
    """Class representing a score for a comparison."""
    genomic_overlap: int = field(default=0)
    matches : int = field(default=0)
    identity: float = field(default=0.0)

    def is_betterEq_than(self, other) -> bool:
        """Compare two scores based on genomic overlap and identity."""
        if other.genomic_overlap == 0:
            return True
        if self.genomic_overlap == 0:
            return False
        if self.identity == other.identity :
                if self.matches == other.matches :
                    return self.genomic_overlap >= other.genomic_overlap
                else:
                    return self.matches >= other.matches
        else:
            return self.identity >= other.identity
    @classmethod
    def add(cls, score1: 'MatchScore', score2: 'MatchScore') -> 'MatchScore':
        """Combine two MatchScore instances."""
        return cls(
            genomic_overlap=score1.genomic_overlap + score2.genomic_overlap,
            identity=score1.identity + score2.identity
        )
    @classmethod
    def max(cls, score1: 'MatchScore', score2: 'MatchScore') -> 'MatchScore':
        """Return the maximum of two MatchScore instances."""
        if score1.is_betterEq_than(score2):
            return score1
        else:
            return score2
    @classmethod
    def max3(cls, score1: 'MatchScore', score2: 'MatchScore', score3: 'MatchScore') -> 'MatchScore':
        """Return the maximum of three MatchScore instances."""
        return cls.max(cls.max(score1, score2), score3)

@define(slots=True, frozen=True, eq=True)
class MismatchInfo:
    """Class representing a mismatch zone between two annotations' loci."""
    zones: list = field(factory=list)  # List of mismatch zone coordinates with included bounds
    nb: int = field(default=0, init=False)  # Number of mismatches in the zone, calculated automatically
    
    def __attrs_post_init__(self):
        """Calculate number of mismatches after initialization."""
        if len(self.zones) % 2 != 0:
            raise ValueError("Mismatch zones must be defined by pairs of coordinates.")
            
        total = 0
        for i in range(0, len(self.zones), 2):
            total += self.zones[i+1] - self.zones[i] + 1
        
        object.__setattr__(self, 'nb', total)
    
    def __str__(self):
        """String representation of the mismatch information."""
        return f"Number of mismatches: {self.nb}, Zones: {self.zones}"
        
    @classmethod
    def create(cls, zone_bounds=None):
        """Factory method to create a MismatchInfo instance."""
        if zone_bounds is None:
            zone_bounds = []
        return cls(zones=zone_bounds.copy())

@define(slots=True, frozen=True, eq=True)
class MrnaMatchInfo:
    """Class representing the result of a comparison between two annotations' loci."""
     # Comparison metrics
    matches: int = field(default=0)
    mismatches_EI: MismatchInfo = field(factory=MismatchInfo)  # Exon-Intron mismatches
    mismatches_RF: MismatchInfo = field(factory=MismatchInfo)  # Exon-Intron mismatches
    genomic_overlap: int = field(default=0)  # Overlap in genomic coordinates
    ref_id: str = field(default='_') 
    alt_id: str =  field(default='_')
    score : MatchScore = field(factory=MatchScore, init=False)  # Score for the match
    
        
    def __attrs_post_init__(self):
        """Calculate identity after initialization."""
        total = self.matches + self.mismatches_EI.nb + self.mismatches_RF.nb
        if total == 0:
            object.__setattr__(self, 'score', MatchScore(self.genomic_overlap, 0, 0))
        else:
            object.__setattr__(self, 'score', MatchScore(self.genomic_overlap, self.matches, self.matches / total))

    def get_identity(self) -> float:
        return self.score.identity

    def has_better_identity_than(self, other: 'MrnaMatchInfo') -> bool:
        return self.score.is_betterEq_than(other.score)
