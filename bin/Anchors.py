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


class AnchorTargetDistances(dict):
    """
    Dictionary object to keep track of anchors, targets, and their distances
    """
    def __init__(self):
        self = {}

    def update_target_distance(self, anchor, target, distance):
        """
        Update {anchor : {target : distance}}
        """
        if anchor not in self:
            self[anchor] = {target : distance}
        else:
            if target not in self[anchor]:
                self[anchor].update({target : distance})
        return self

    def update_topTargets_distances(self, anchor, topTargets_distances):
        """
        Update {anchor : {target : distance}} with the distances of the topTargets
        """
        if anchor in self:
            self[anchor].update(dict(topTargets_distances))
        else:
            self[anchor] = dict(topTargets_distances)
        return self

    def has_distance(self, anchor, target):
        """
        Return true if (anchor, target) has a pre-computed distance
        """
        if target in self[anchor]:
            return True
        else:
            return False

    def get_distance(self, anchor, target):
        """
        Return distance for anchor:target pair
        """
        return self[anchor][target]


class AnchorTargetsSamples(dict):
    def __init__(self):
        self = {}

    def update(self, anchor, target, sample):
        """
        Update {anchor : {sample : {target : count}}
        """
        if anchor in self:
            if sample in self[anchor]:
                if target in self[anchor][sample]:
                    self[anchor][sample][target] += 1
                else:
                    self[anchor][sample][target] = 1
            else:
                self[anchor][sample] = {target : 1}
        else:
            self[anchor] = {sample : {target : 1}}

        return self

    def get_counts_df(self, final_anchors_list):
        """
        Return a table of anchors and their target counts
        """
        # nested dict to dataframe
        df_rows = []
        for anchor, anchor_dict in self.items():
            for sample, target_dict in anchor_dict.items():
                for target, count in target_dict.items():
                    df_rows.append([anchor, sample, target, count])

        # pivot table to anchors x counts
        df = (
            pd.DataFrame(
                df_rows,
                columns=['anchor', 'sample', 'target', 'count']
            )
            .pivot(
                index=['anchor', 'target'],
                columns=['sample'],
                values='count'
            )
            .reset_index()
        )
        df = df[df['anchor'].isin(final_anchors_list)]
        return df


    def get_anchor_counts_df(self, anchor):
        """
        Return a table of anchors x samples counts for a specific anchor, sorted by abundance
        """
        # filter dict for the anchor
        d = {k:v for (k,v) in self.items() if k == anchor}

        # nested dict to table
        df_rows = []
        for anchor, anchor_dict in d.items():
            for sample, target_dict in anchor_dict.items():
                for target, count in target_dict.items():
                    df_rows.append([anchor, sample, target, count])

        # pivot table to anchors x counts
        df = (
            pd.DataFrame(
                df_rows,
                columns=['anchor', 'sample', 'target', 'count']
            )
            .pivot(
                index=['anchor', 'target'],
                columns=['sample'],
                values='count'
            )
        )
        df.columns.name = None

        # remove anchor column and set targets as index
        anchor_df = (
            df
            .reset_index()
            .drop('anchor', axis=1)
            .set_index('target')
        )

        # sort targets by abundance
        anchor_df['count'] = anchor_df.sum(axis=1)
        anchor_df = (
            anchor_df
            .sort_values(
                'count',
                ascending=False)
            .drop(
                'count',
                axis=1
            )
            .reset_index()
        )

        return anchor_df


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
    def __init__(self, anchor_counts, anchor_targets_samples, anchor_targets, anchor_scores_topTargets, anchor_target_distances, anchor_status):
        self.ignorelist = {}
        self.keeplist = {}
        self.anchor_counts = anchor_counts
        self.anchor_targets_samples = anchor_targets_samples
        self.anchor_targets = anchor_targets
        self.anchor_scores_topTargets = anchor_scores_topTargets
        self.anchor_target_distances = anchor_target_distances
        self.anchor_status = anchor_status

    def update_ignorelist(self, anchor, read_counter_freeze):
        """
        Updates ignorelist and removes anchor from all of our dictionaries
        """
        # only add anchor to ignorlist if not in read_counter_freeze or if ignorelist is not larger than 4M
        if (not read_counter_freeze) and (len(self.ignorelist) < 4000000):
            if anchor not in self.ignorelist:
                self.ignorelist[anchor] = True

        # pop anchor from all dicts
        self.anchor_counts.total_counts.pop(anchor, None)
        self.anchor_counts.all_target_counts.pop(anchor, None)
        self.anchor_counts.top_target_counts.pop(anchor, None)
        self.anchor_targets_samples.pop(anchor, None)
        self.anchor_targets.pop(anchor, None)
        self.anchor_scores_topTargets.pop(anchor, None)
        self.anchor_target_distances.pop(anchor, None)
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



class AnchorScoresTopTargets(dict):
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

    def compute_target_distance(self, anchor, target):
        """
        Return the min distance of a target to any of its anchor's topTargets
        """
        topTargets = self[anchor].index.tolist()
        distance = min([utils.get_distance(target, t) for t in topTargets])
        return distance

    def get_summary_scores(self, group_ids_dict, use_std):
        """
        Return scores, scaled and sorted by abundance
        """
        # Convert dict to dataframe
        scores_df = pd.DataFrame()
        for anchor, (scores, _) in self.items():
            scores_df = scores_df.append(pd.Series(scores, name=anchor))

        # drop any anchors that do not have scores for at least 3 samples
        scores_df = scores_df.dropna(thresh=3)

        # if we are using standard deviation score or there is no metadata component, get the std of scores
        if (use_std or len(set(group_ids_dict.values())) == 1):
            df = (
                pd.DataFrame(
                    scores_df.std(axis=1),
                    columns=['score']
                )
            )

        # else, compute the anchor significane score as sum(C_i*S_i) - sum(S_i(sum(C_i) / n_c))
        else:

            # sum(C_i*S_i)
            first_sum = (
                scores_df
                .assign(**group_ids_dict)
                .mul(scores_df)
                .dropna(axis=1)
                .sum(axis=1)
            )

            # sum(S_i(sum(C_i) / n_c))
            second_sum = (
                scores_df
                .apply(
                    lambda row:
                    (
                        (
                            sum(
                                ([group_ids_dict[col] for col in scores_df.columns if not np.isnan(row[col])])
                            )
                            *
                            sum(
                                [row[col] for col in scores_df.columns if not np.isnan(row[col])]
                            )
                        )
                        /
                        len(scores_df.columns)
                    ),
                    axis=1
                )
            )

            df = (
                pd.DataFrame(
                    first_sum - second_sum,
                    columns=['score']
                )
                .sort_values(
                    'score',
                    ascending=False
                )
            )
        return df

