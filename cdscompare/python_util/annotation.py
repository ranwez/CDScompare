# -*- coding: utf-8 -*-

from pathlib import Path
from typing import List
from attrs import define


@define(frozen=True)
class Annotation:
    id: str
    path: Path


@define(frozen=True)
class AnnotationPair:
    ref: Annotation
    alt: Annotation

    def output_filenames(self, out_dir: Path) -> tuple[Path, Path]:
        base = f"{self.ref.id}_vs_{self.alt.id}"
        return (
            out_dir / f"{base}.csv",
            out_dir / f"{base}.txt",
        )


@define
class AnnotationSet:
    ref: Annotation
    alts: List[Annotation]

    @classmethod
    def from_paths(cls, ref_path: Path, alt_paths: List[Path]) -> "AnnotationSet":
        ref = cls._make_annotation(ref_path)
        alts = [cls._make_annotation(p) for p in alt_paths]
        instance = cls(ref=ref, alts=alts)
        instance._check_unique_ids()
        return instance

    def pairs(self) -> list[AnnotationPair]:
        return [AnnotationPair(self.ref, alt) for alt in self.alts]

    def synthesis_filename(self, out_dir: Path) -> Path:
        return out_dir / f"synthesis_{self.ref.id}.csv"

    @staticmethod
    def _make_annotation(path: Path) -> Annotation:
        return Annotation(id=path.stem, path=path)

    def _check_unique_ids(self) -> None:
        ids = [self.ref.id] + [a.id for a in self.alts]
        duplicates = {i for i in ids if ids.count(i) > 1}
        if duplicates:
            raise ValueError(
                f"Duplicate annotation identifiers detected: {', '.join(sorted(duplicates))}. "
                "Please rename the input GFF files so that each has a unique base name."
            )
