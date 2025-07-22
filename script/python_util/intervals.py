#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: vetea, ranwez
"""
import array
from itertools import chain
## This class represents a list of sorted, non overlapping intervals and implements multiple basic
# operations between two such lists (intersection, union, difference, symmetric difference).
# [[start1, end1] [start2, end2] ...] is represented as an array of integers [start1, end1, start2, end2, ...]
#
# @remark The logic is strongly inspired from https://stackoverflow.com/a/20062829 
class OrderedIntervals:
    
    def __init__(self, intervals:array, include_ub=False):
        if include_ub:
            self.intervals = OrderedIntervals.transform_intervals_to_exclude_ub(intervals)
        else:
            self.intervals = intervals if isinstance(intervals, array.array) else array.array('L', intervals)
    
    @staticmethod
    def transform_intervals_to_exclude_ub(intervals: array) -> array.array:
        """Transform [start,end] intervals into equivalent [start, end+1["""
        transformed_intervals = array.array('L', [0] * len(intervals))
        for i in range(0, len(intervals), 2):
            transformed_intervals[i] = intervals[i]
            transformed_intervals[i+1] = intervals[i + 1] + 1
        return transformed_intervals
    
    def as_list_with_included_ub(self) -> list[int]:
        """Transform [start,end+1[ intervals into equivalent [start, end] intervals"""
        return [val if i % 2 == 0 else val - 1 
                for i, val in enumerate(self.intervals)]

    def total_length(self) -> int:
        return sum(self.intervals[i + 1] - self.intervals[i] for i in range(0, len(self.intervals), 2))

    @staticmethod
    def new(intervals, include_ub=False):
        return OrderedIntervals(intervals,include_ub)
    
    def union(self, other):
        return self.merge(other, lambda a, b: a or b)

    def inter_union_symdiff(self, other):
        return self.triple_merge(other, lambda a, b: a and b, lambda a, b: a or b, lambda a, b: a ^ b)

    def intersection(self, other):
        return self.merge(other, lambda a, b: a and b)
    
    def difference(self, other):
        """Returns the intervals of the first list not present in the second"""
        return self.merge(other, lambda a, b: a and not b)

    def symmetric_difference(self, other):
        return self.merge(other, lambda a, b: a ^ b)
    
    def merge(self, other, keep_operator) -> 'OrderedIntervals':
        """ The key method used to merge two OrderedIntervals objects."""
        if not self.intervals and not other.intervals:
            return OrderedIntervals(array.array('L'))

        max_size = len(self.intervals) + len(other.intervals)
        res = array.array('L', [0] * max_size)  
        res_idx = 0  

        sentinel = max(self.intervals[-1] if self.intervals else 0, 
                     other.intervals[-1] if other.intervals else 0) + 1

        iter0 = chain(self.intervals, [sentinel])
        iter1 = chain(other.intervals, [sentinel])
        bound0 = next(iter0)
        bound1 = next(iter1)
        is_lb0 = True
        is_lb1 = True

        scan = min(bound0, bound1)
        next_res_is_lb = True

        while scan < sentinel:
            in0 = (scan >= bound0) == is_lb0
            in1 = (scan >= bound1) == is_lb1
            in_res = keep_operator(in0, in1)

            if in_res == next_res_is_lb:
                res[res_idx] = scan
                res_idx += 1
                next_res_is_lb = not next_res_is_lb

            if scan == bound0:
                bound0 = next(iter0)
                is_lb0 = not is_lb0
            if scan == bound1:
                bound1 = next(iter1)
                is_lb1 = not is_lb1

            scan = min(bound0, bound1)

        final_res = array.array('L', res[:res_idx])
        return OrderedIntervals(final_res)

    def triple_merge(self, other, op1, op2, op3) -> tuple['OrderedIntervals', 'OrderedIntervals', 'OrderedIntervals']:
        """Merge two OrderedIntervals objects with three different operations."""
        if not self.intervals and not other.intervals:
            empty = OrderedIntervals(array.array('L'))
            return empty, empty, empty

        max_size = len(self.intervals) + len(other.intervals)
        
        res1 = array.array('L', [0] * max_size)
        res2 = array.array('L', [0] * max_size)
        res3 = array.array('L', [0] * max_size)
        
        idx1, idx2, idx3 = 0, 0, 0

        sentinel = max(self.intervals[-1] if self.intervals else 0, 
                      other.intervals[-1] if other.intervals else 0) + 1

        iter0 = chain(self.intervals, [sentinel])
        iter1 = chain(other.intervals, [sentinel])
        bound0 = next(iter0)
        bound1 = next(iter1)
        is_lb0 = True
        is_lb1 = True

        scan = min(bound0, bound1)
        next_is_lb1 = True
        next_is_lb2 = True
        next_is_lb3 = True

        while scan < sentinel:
            in0 = (scan >= bound0) == is_lb0
            in1 = (scan >= bound1) == is_lb1

            if op1(in0, in1) == next_is_lb1:
                res1[idx1] = scan
                idx1 += 1
                next_is_lb1 = not next_is_lb1
            
            if op2(in0, in1) == next_is_lb2:
                res2[idx2] = scan
                idx2 += 1
                next_is_lb2 = not next_is_lb2
            
            if op3(in0, in1) == next_is_lb3:
                res3[idx3] = scan
                idx3 += 1
                next_is_lb3 = not next_is_lb3

            if scan == bound0:
                bound0 = next(iter0)
                is_lb0 = not is_lb0
            if scan == bound1:
                bound1 = next(iter1)
                is_lb1 = not is_lb1

            scan = min(bound0, bound1)

        final_res1 = array.array('L', res1[:idx1])
        final_res2 = array.array('L', res2[:idx2])
        final_res3 = array.array('L', res3[:idx3])

        return OrderedIntervals(final_res1), OrderedIntervals(final_res2), OrderedIntervals(final_res3)
        