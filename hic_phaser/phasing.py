import logging
from dataclasses import dataclass
from typing import Dict
from .gfa import GfaSegments


@dataclass
class LinkCounts:
    counts: Dict[tuple, int]  # counts on "links" between segments


class HicPhaser:
    def __init__(self, link_counts: LinkCounts):
        pass


