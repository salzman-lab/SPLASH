#!/usr/bin/env python3

from Bio import SeqIO


class Read():
    """
    Object for a fastq read
    """
    def __init__(self, read):
        self.read = read

    def get_anchors(self, anchor_mode, window_slide, lookahead, kmer_size):
        """
        Get list of chunked anchors from read
        """
        if anchor_mode == 'chunk':
            step_size = kmer_size
        elif anchor_mode == 'tile':
            step_size = window_slide

        last_base = len(self.read) - (lookahead + 2 * kmer_size)
        anchor_list = [
            self.read[0+i:kmer_size+i]
            for i
            in range(0, last_base, step_size)
            if len(self.read[0+i:kmer_size+i])==kmer_size
        ]
        return anchor_list

    def get_target(self, anchor, lookahead, kmer_size):
        """
        Get targets for a given anchor that are size=target_len and target_dist away from the anchor
        """
        # get anchor position
        anchor_end = self.read.index(anchor) + len(anchor)

        # get target position
        target_start = anchor_end + lookahead
        target_end = target_start + kmer_size

        target = self.read[target_start:target_end]

        # if adj anchor exists, add adj anchor to anchor_dict
        if len(target) == kmer_size:
            return target
        else:
            return None
