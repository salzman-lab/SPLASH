#!/usr/bin/env python

from Bio import SeqIO


class Read():
    """
    Object for a fastq read
    """
    def __init__(self, read):
        self.read = read
    
    def get_anchors(self, anchor_len, anchor_mode, window_slide):
        """
        Get list of chunked anchors from read
        """
        if anchor_mode == 'chunk':
            step_size = anchor_len 
        elif anchor_mode == 'tile':
            step_size = window_slide
        anchor_list = [
            self.read[0+i:anchor_len+i] 
            for i 
            in range(0, len(self.read), step_size) 
            if len(self.read[0+i:anchor_len+i])==anchor_len
        ]
        return anchor_list

    def get_target(self, anchor, target_dist, target_len):
        """
        Get targets for a given anchor that are size=target_len and target_dist away from the anchor
        """
        # get anchor position
        anchor_end = self.read.index(anchor) + len(anchor)

        # get target position
        target_start = anchor_end + target_dist
        target_end = target_start + target_len

        target = self.read[target_start:target_end]

        # if adj anchor exists, add adj anchor to anchor_dict
        if len(target) == target_len:
            return target
        else:
            return None