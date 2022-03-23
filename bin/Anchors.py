#!/usr/bin/env python

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
    List object that contains anchors if they are in phase 2
    """
    def __init__(self):
        self.anchors = []

    def assign_phase_2(self, anchor):
        self.anchors.append(anchor)
        return self

    def is_phase_2(self, anchor):
        if anchor in self.anchors:
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
        # nested dict to table
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

        anchor_df = df.reset_index().drop('anchor', axis=1).set_index('target')
        anchor_df['count'] = anchor_df.sum(axis=1)
        anchor_df = anchor_df.sort_values('count', ascending=False).drop('count', axis=1).reset_index()
        return anchor_df

    def reached_phase_2_threshold(self, anchor, target, sample):
        """
        Returns true if (anchor:target:sample) has at least 10 phase_2 scores
        """
        if self[anchor][sample][target] >= 20:
            return True
        else:
            return False


class AnchorCounts(dict):
    """
    Dictionary object of anchors and their total counts
    """
    def __init__(self, num_samples):
        self.counter = {}
        self.num_samples = num_samples
        self.num_unique_anchors = 0

    def contains(self, anchor):
        """
        Returns true if anchor is in dict
        """
        if anchor in self.counter:
            return True
        else:
            return False

    def update(self, anchor, iteration):
        """
        Update {anchor : count}
        """

        if anchor not in self.counter:
            # if this is a new anchor, increment relative to sample number and iteration
            L = math.floor(iteration / 100000)
            increment = self.num_samples * L
            self.counter[anchor] = 1 + increment

            # update number of unique anchors count
            self.num_unique_anchors += 1

        else:
            # if not new anchor, increment relative to sample number and iteration
            increment = self.num_samples
            self.counter[anchor] += increment


    def get_count(self, anchor):
        """
        Get count for an anchor
        """
        return self.counter[anchor]

    def get_ignorelist_anchors(self, anchor_min_count):
        """
        Return anchors that do not have at least anchor_min_count number of counts
        """
        return [anchor for anchor, count in self.counter.items() if count < anchor_min_count]


class AnchorScoresTopTargets(dict):
    """
    Dictionary object to keep track of anchors, their current scores, and their top targets
    """
    def __init__(self):
        self = {}

    def update(self, anchor, scores, top_targets):
        """
        Update {anchor : [scores, top_targets]}
        """
        self[anchor] = [scores, top_targets]
        return self

    def get_score(self, anchor):
        """
        Return score for a given anchor
        """
        return self[anchor][0]

    def get_topTargets(self, anchor):
        """
        Return topTargets for a given anchor
        """
        return self[anchor][1]

    def compute_target_distance(self, anchor, target):
        """
        Return the min distance of a target to any of its anchor's topTargets
        """
        topTargets = self.get_topTargets(anchor)
        distance = min([utils.get_distance(target, t) for t in topTargets])
        return distance

    def get_target_scores(self, final_anchors_list):
        """
        Return a df of anchors, targets, and their current scores
        """
        scores_df = pd.DataFrame()
        for anchor, (scores, _) in self.items():
            scores_df = scores_df.append(pd.Series(scores, name=anchor))
        scores_df.index.name = 'anchor'
        scores_df = scores_df.reset_index()
        scores_df = scores_df[scores_df['anchor'].isin(final_anchors_list)]
        return scores_df

    def get_phase_scores_df(self):
        """Temp helper func"""
        scores_df = pd.DataFrame()
        for anchor, (scores, _) in self.items():
            scores_df = scores_df.append(pd.Series(scores, name=anchor))
        scores_df.index.name = 'anchor'

        return scores_df.reset_index()

    def get_num_scores(self):
        """
        Return the number of anchors with scores
        """
        scores_df = pd.DataFrame()
        for anchor, (scores, _) in self.items():
            scores_df = scores_df.append(pd.Series(scores, name=anchor))

        scores_df = scores_df.dropna(thresh=3)
        return len(scores_df)

    def get_summary_scores(self, group_ids_dict, use_std):
        """
        Return scores, scaled and sorted by abundance
        """
        scores_df = pd.DataFrame()
        for anchor, (scores, _) in self.items():
            scores_df = scores_df.append(pd.Series(scores, name=anchor))

        scores_df = scores_df.dropna(thresh=3)

        if (use_std or len(set(group_ids_dict.values())) == 1):
            df = (
                pd.DataFrame(
                    scores_df.std(axis=1),
                    columns=['score']
                )
            )
        else:
            first_sum = (
                scores_df
                .assign(**group_ids_dict)
                .mul(scores_df)
                .dropna(axis=1)
                .sum(axis=1)
            )

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

    def get_targets(self, anchor):
        """
        Return the unique targets of an anchor
        """
        return self[anchor]


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

    def update_ignorelist(self, anchor, read_counter_freeze, anchor_counter_freeze):
        """
        Updates ignorelist and removes anchor from all of our dictionaries
        """
        if not any([read_counter_freeze, anchor_counter_freeze]):
            if anchor not in self.ignorelist:
                self.ignorelist[anchor] = True

        self.anchor_counts.pop(anchor, None)
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
