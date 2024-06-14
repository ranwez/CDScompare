class OrderedIntervals:
    def __init__(self, intervals, include_ub=False):
        if include_ub:
            self.intervals = OrderedIntervals.transform_intervals_to_exclude_ub(intervals)
        else:
            self.intervals = intervals
    
    @staticmethod
    def transform_intervals_to_exclude_ub(intervals):
        transformed_intervals = []
        for i in range(0, len(intervals), 2):
            lower_bound = intervals[i]
            upper_bound = intervals[i + 1] + 1
            transformed_intervals.extend([lower_bound, upper_bound])
        return(transformed_intervals)

   
    def get_intervals_with_included_ub(self):
        transformed_intervals = []
        for i in range(0, len(self.intervals), 2):
            lower_bound = self.intervals[i]
            upper_bound = self.intervals[i + 1] - 1
            transformed_intervals.extend([lower_bound, upper_bound])
        return(transformed_intervals)


    def total_length(self):
        return sum(self.intervals[i + 1] - self.intervals[i] for i in range(0, len(self.intervals), 2))

    @staticmethod
    def new(intervals,include_ub=False):
        return OrderedIntervals(intervals,include_ub)

    def union(self, other):
        res =self.merge(other, lambda a, b: a or b).get_intervals_with_included_ub()

    def union(self, other):
        return self.merge(other, lambda a, b: a or b)

    def intersection(self, other):
        return self.merge(other, lambda a, b: a and b)

    def difference(self, other):
        return self.merge(other, lambda a, b: a and not b)

    def symmetric_difference(self, other):
        return self.merge(other, lambda a, b: a ^ b)

    def merge(self, other, keep_operator):
        res = []
        if not self.intervals and not other.intervals:
            return OrderedIntervals(res)

        sentinel = max(self.intervals[-1], other.intervals[-1]) + 1

        interval_iters = (
            iter(self.intervals + [sentinel]),
            iter(other.intervals + [sentinel]),
        )

        current_bound = (
            next(interval_iters[0]),
            next(interval_iters[1]),
        )
        scan = min(current_bound[0], current_bound[1])
        current_is_lb = (True, True)
        next_res_is_lb = True

        while scan < sentinel:
            in_interval = (
                (scan >= current_bound[0]) == current_is_lb[0],
                (scan >= current_bound[1]) == current_is_lb[1],
            )
            in_res = keep_operator(in_interval[0], in_interval[1])

            if in_res == next_res_is_lb:
                res.append(scan)
                next_res_is_lb = not next_res_is_lb

            if scan == current_bound[0]:
                current_bound = (next(interval_iters[0]), current_bound[1])
                current_is_lb = (not current_is_lb[0], current_is_lb[1])

            if scan == current_bound[1]:
                current_bound = (current_bound[0], next(interval_iters[1]))
                current_is_lb = (current_is_lb[0], not current_is_lb[1])

            scan = min(current_bound[0], current_bound[1])

        return OrderedIntervals(res)


# Example usage
intervals1 = OrderedIntervals.new([1, 2, 8, 10, 16, 18])
intervals2 = OrderedIntervals.new([2, 4, 7, 10, 15, 16, 18, 20])
intervals3 = intervals1.union(intervals2)
print(intervals3.intervals)  # Output: [1, 4, 7, 10, 15, 20]
assert intervals3.intervals == [1, 4, 7, 10, 15, 20]
lg = intervals1.total_length()
print(lg)  # Output: 5

intervals1bis = OrderedIntervals.new([1, 2, 8, 10, 16, 18], include_ub=True)
lg1bis = intervals1bis.total_length()
print(lg1bis)  # Output: 7

intervals4 = intervals1.intersection(intervals2)
print(intervals4.intervals)  # Output: [2, 10, 16, 18]
assert intervals4.intervals == [8, 10]
