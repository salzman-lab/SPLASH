#!/usr/bin/env python3

import pandas as pd
import gzip
import os
import re
import sys
import argparse
import logging
import numpy as np
from Bio import SeqIO
import nltk
import math
import logging
import utils


class AnchorStatus(dict):
    """
    Dictionary object that contains anchors if they are in phase 2
    """
    def __init__(self):
        self = {}

    def assign_phase_2(self, anchor):
        """
        Add {anchor: True} if an anchor is assigned to phase_2
        """
        self[anchor] = True
        return self

    def is_phase_2(self, anchor):
        """
        Returns True if an anchor is in phase_2
        """
        if anchor in self:
            return True
        else:
            return False


class AnchorTopTargetsDistances(dict):
    """
    Dictionary object to keep track of anchors, targets, and their distances
    """
    def __init__(self):
        self.distances = {}
        self.mu = {}

    def update_distances(self, anchor, topTargets_distances):
        """
        Update {anchor : {target : distance}} with the distances of the topTargets
        """
        if anchor in self.distances:
            self.distances[anchor].update(dict(topTargets_distances))
        else:
            self.distances[anchor] = dict(topTargets_distances)
        return self.distances

    def update_mu(self, anchor, mu):
        """
        Updte {anchor: mu}
        """
        self.mu[anchor] = mu
        return self.mu

    def has_distance(self, anchor, target):
        """
        Return true if (anchor, target) has a pre-computed distance
        """
        if target in self.distances[anchor]:
            return True
        else:
            return False

    def get_distance(self, anchor, target):
        """
        Return distance for anchor:target pair
        """
        return self.distances[anchor][target]

    def get_mu(self, anchor):
        """
        Return distance for anchor:target pair
        """
        return self.mu[anchor]

class AnchorTargetsSamples():
    def __init__(self, sample_index_dict):
        self.counter = {}
        self.sample_index_dict = sample_index_dict
        self.name_dict = {v:k for k,v in self.sample_index_dict.items()}

    def update(self, anchor, target, sample):
        """
        Update {anchor : {target : [sample_1_count, sample_2_count, ...]}
        """
        # get sample index
        sample_index = self.sample_index_dict[sample]

        # intialise
        if anchor not in self.counter:
            self.counter[anchor] = {target : [0] * len(self.sample_index_dict)}
            self.counter[anchor][target][sample_index] = 1

        else:
            if target not in self.counter[anchor]:
                self.counter[anchor].update({target : [0] * len(self.sample_index_dict)})
                self.counter[anchor][target][sample_index] = 1
            else:
                self.counter[anchor][target][sample_index] += 1

    def get_anchor_counts_df(self, anchor):
        """
        Return a table of anchors x samples counts for a specific anchor, sorted by abundance
        """
        df = pd.DataFrame(self.counter[anchor]).T
        df.columns = [self.name_dict[c] for c in df.columns]
        df['count'] = df.sum(axis=1)
        df = (
            df
            .sort_values('count', ascending=False)
            .drop('count', axis=1)
            .reset_index()
            .rename(columns={'index':'target'})
        )
        return df


class AnchorCounts(dict):
    """
    Dictionary object of anchors and their total counts
    """
    def __init__(self, num_samples):
        self.total_counts = {}
        self.all_target_counts = {}
        self.top_target_counts = {}
        self.num_samples = num_samples
        self.num_unique_anchors = 0

    def contains(self, anchor):
        """
        Returns true if anchor is in dict
        """
        if anchor in self.total_counts:
            return True
        else:
            return False

    def update_total_counts(self, anchor, iteration):
        """
        Update {anchor : count}
        """
        if anchor not in self.total_counts:
            # if this is a new anchor, increment relative to sample number and iteration
            L = math.floor(iteration / 100000)
            increment = self.num_samples * L
            self.total_counts[anchor] = 1 + increment

            # update number of unique anchors count
            self.num_unique_anchors += 1

        else:
            # if not new anchor, increment relative to sample number and iteration
            increment = self.num_samples
            self.total_counts[anchor] += increment

    def update_all_target_counts(self, anchor, sample, iteration):
        """
        Update {anchor : {sample : count}} relative to the iteration and num_samples
        """
        L = math.floor(iteration / 100000)

        init_count = self.num_samples
        increment_count = 1 + (self.num_samples * L)

        # if the anchor is already in the dict
        if anchor in self.all_target_counts:
            if sample in self.all_target_counts[anchor]:
                # if this is not new, increment relative to sample number and iteration
                self.all_target_counts[anchor][sample] += increment_count

            else:
                # if this is new, intialize with init_count
                self.all_target_counts[anchor][sample] = init_count

        # if this is a new anchor, initialize with int_count
        else:
            self.all_target_counts[anchor] = {sample : init_count}

    def get_total_counts(self, anchor):
        """
        Get count for an anchor
        """
        return self.total_counts[anchor]

    def get_all_target_counts(self, anchor, sample):
        """
        Get count for an anchor
        """
        return self.all_target_counts[anchor][sample]

    def get_ignorelist_anchors(self, anchor_min_count):
        """
        Return anchors that do not have at least anchor_min_count number of counts
        """
        return [anchor for anchor, count in self.total_counts.items() if count < anchor_min_count]


class AnchorTargets(dict):
    """
    Dictionary object to keep track of anchors and their targets
    """
    def __init__(self):
        self = {}

    def update_target(self, anchor, target):
        """
        Add unique target to {anchor : target}
        """
        if anchor in self:
            # if this is a new target, add target to anchor
            if target not in self[anchor]:
                self[anchor].append(target)
        else:
            self[anchor] = [target]
        return self

    def num_targets(self, anchor):
        """
        Return the number of unique targets encountered for that anchor
        """
        return len(self[anchor])

    def is_topTarget(self, anchor, target):
        """
        Return the unique targets of an anchor
        """
        if anchor in self:
            if target in self[anchor]:
                return True
            else:
                return False
        else:
            return False


class StatusChecker():
    """
    Object to keep track of keeplists and ignorelists
    """
    def __init__(self, anchor_counts, anchor_targets_samples, anchor_targets, anchor_topTargets_scores, anchor_topTargets_distances, anchor_status, max_ignorelist):
        self.ignorelist = {}
        self.keeplist = {}
        self.anchor_counts = anchor_counts
        self.anchor_targets_samples = anchor_targets_samples
        self.anchor_targets = anchor_targets
        self.anchor_topTargets_scores = anchor_topTargets_scores
        self.anchor_topTargets_distances = anchor_topTargets_distances
        self.anchor_status = anchor_status
        self.max_ignorelist = max_ignorelist

    def update_ignorelist(self, anchor, read_counter_freeze):
        """
        Updates ignorelist and removes anchor from all of our dictionaries
        """
        # only add anchor to ignorlist if not in read_counter_freeze or if ignorelist is not larger than 4M
        if (not read_counter_freeze) and (len(self.ignorelist) < self.max_ignorelist):
            if anchor not in self.ignorelist:
                self.ignorelist[anchor] = True

        # pop anchor from all dicts
        self.anchor_counts.total_counts.pop(anchor, None)
        self.anchor_counts.all_target_counts.pop(anchor, None)
        self.anchor_counts.top_target_counts.pop(anchor, None)
        self.anchor_targets_samples.counter.pop(anchor, None)
        self.anchor_targets.pop(anchor, None)
        self.anchor_topTargets_scores.pop(anchor, None)
        self.anchor_topTargets_distances.distances.pop(anchor, None)
        self.anchor_topTargets_distances.mu.pop(anchor, None)
        self.anchor_status.pop(anchor, None)

    def update_keeplist(self, anchor):
        """
        Updates keeplist
        """
        if anchor not in self.keeplist:
            self.keeplist[anchor] = True


    def is_ignorelisted(self, anchor):
        """
        Returns true if an anchor has previously ignorelisted
        """
        if anchor in self.ignorelist:
            return True
        else:
            return False



class AnchorTopTargetsScores(dict):
    """
    Dictionary object to keep track of anchors, their current scores, and their top targets
    """
    def __init__(self):
        self = {}

    def initialise(self, anchor, scores):
        """
        Add {anchor : sample_scores}
        """
        self[anchor] = scores
        return self

    def update(self, anchor, sample, score):
        """
        Update {anchor : sample_scores}
        """
        self[anchor][sample] = score
        return self

    def get_score(self, anchor, sample):
        """
        Return score for a given anchor
        """
        if sample in self[anchor]:
            return self[anchor][sample]
        else:
            return 0

    def compute_target_distance(self, anchor, target, distance_type):
        """
        Return the min distance of a target to any of its anchor's topTargets
        """
        topTargets = self[anchor].index.tolist()
        distance = min([utils.get_distance(target, t, distance_type) for t in topTargets])
        return distance

    def get_blacklist_anchors(self, use_std):
        """
        Get sum of all scores and sort by sum descending
        """
        df = pd.DataFrame(self).T

        if use_std:
            df['z'] = (
                df
                .std(axis=1)
                .abs()
            )
        else:
            df['z'] = (
                df
                .sum(axis=1)
                .abs()
            )

        anchors = (
            df[df['z'] < df['z'].quantile(0.2)]
            .index
            .to_list()
        )
        return anchors

    def get_final_anchors(self, num_keep_anchors, use_std):
        """
        Get final anchors for parse_anchors
        """
        if len(self) < num_keep_anchors:
            anchors = list(self.keys())
        else:
            df = pd.DataFrame(self).T
            if use_std:
                df['z'] = (
                    df
                    .std(axis=1)
                    .abs()
                )
            else:
                df['z'] = (
                    df
                    .sum(axis=1)
                    .abs()
                )

            anchors = (
                df['z']
                .sort_values(ascending=False)
                .head(num_keep_anchors)
                .index
                .to_list()
            )
        return anchors

