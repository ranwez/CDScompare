# -*- coding: utf-8 -*-

import cdscompare.python_util.intervals as iu


# test function for the 'intervals.py' class method 'transform_intervals_to_exclude_ub'
def test_transform_intervals_to_exclude_ub():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        "basic" : [[100, 129, 150, 209, 240, 299], 
                   [100, 130, 150, 210, 240, 300]],
        
        "identical" : [[100, 129, 150, 209, 240, 299], 
                   [100, 130, 150, 210, 240, 300]],
        
        "minus-CDS" : [[100, 129, 240, 299], 
                   [100, 130, 240, 300]],
        
        "fusion" : [[100, 209, 240, 299], 
                   [100, 210, 240, 300]],
        
        "shift" : [[100, 129, 151, 209, 240, 299], 
                   [100, 130, 151, 210, 240, 300]],
        
        "reverse": [[0, 59, 90, 149, 170, 199],
                    [0, 60, 90, 150, 170, 200]],
        
        "diff-start-before" : [[40, 69, 90, 149, 180, 239], 
                               [40, 70, 90, 150, 180, 240]],
        
        "diff-start-after" : [[160, 189, 210, 269, 300, 359], 
                              [160, 190, 210, 270, 300, 360]],
        
        "basic-2-loci (second locus)" : [[600, 699, 800, 899], 
                                         [600, 700, 800, 900]],
        
        "length_computation_ref" : [[8, 13], 
                   [8, 14]],
        
        "length_computation_alt" : [[1, 4, 12, 13], 
                   [1, 5, 12, 14]]
        
        }
    
    print("\n*************Testing the transform_intervals_to_exclude_ub function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        # transform_intervals_to_exclude_ub is called by the class constructor
        interval = iu.OrderedIntervals(test_dict[test][0], True)
        result = interval.intervals
        print(f"result : {result}\n")
        print(test_dict[test][1])
        assert list(result) == test_dict[test][1]
        
        
# test function for the 'intervals_utils.py' class method 'total_length'
def test_total_length():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        "basic" : [[100, 129, 150, 209, 240, 299],
                   150],
        
        "identical" : [[100, 129, 150, 209, 240, 299],
                   150],
        
        "minus-CDS" : [[100, 129, 240, 299],
                   90],
        
        "fusion" : [[100, 209, 240, 299],
                   170],
        
        "shift" : [[100, 129, 151, 209, 240, 299],
                   149],
        
        "reverse" : [[0, 59, 90, 149, 170, 199],
                     150],
        
        "diff-start-before" : [[40, 69, 90, 149, 180, 239],
                               150],
        
        "diff-start-after" : [[160, 189, 210, 269, 300, 359],
                              150],
        
        "basic-2-loci (second locus)" : [[600, 699, 800, 899],
                                         200],
        
        "length_computation_ref" : [[8, 13],
                                    6],
        
        "length_computation_alt" : [[1, 4, 12, 13],
                                    6],
        
        }
    
    print("\n*************Testing the total_length function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        interval = iu.OrderedIntervals(test_dict[test][0], True)
        result = interval.total_length()
        print(f"result : {result}\n")
        print(test_dict[test][1])
        assert result == test_dict[test][1]
        
        
# test function for the 'intervals_utils.py' class method 'intersection'
def test_intersection():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        
        "identical" : [[100, 129, 150, 209, 240, 299],
                       [100, 129, 150, 209, 240, 299],
                   [100, 129, 150, 209, 240, 299]],
        
        "minus-CDS" : [[100, 129, 150, 209, 240, 299],
                       [100, 129, 240, 299],
                   [100, 129, 240, 299]],
        
        "fusion" : [[100, 129, 150, 209, 240, 299],
                    [100, 209, 240, 299],
                   [100, 129, 150, 209, 240, 299]],
        
        "shift" : [[100, 129, 150, 209, 240, 299],
                   [100, 129, 151, 209, 240, 299],
                   [100, 129, 151, 209, 240, 299]],
        
        "reverse" : [[0, 59, 90, 149, 170, 199],
                     [0, 59, 90, 149, 170, 199],
                     [0, 59, 90, 149, 170, 199]],
        
        "reverse-modified" : [[0, 59, 90, 149, 170, 199],
                              [0, 59, 80, 149, 170, 199],
                              [0, 59, 90, 149, 170, 199]],
        
        "diff-start-before" : [[100, 129, 150, 209, 240, 299],
                               [40, 69, 90, 149, 180, 239],
                               [100, 129, 180, 209]],
        
        "diff-start-after" : [[100, 129, 150, 209, 240, 299],
                              [160, 189, 210, 269, 300, 359],
                              [160, 189, 240, 269]],
        
        }
    
    print("\n*************Testing the intersection function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        interval1 = iu.OrderedIntervals(test_dict[test][0], True)
        interval2 = iu.OrderedIntervals(test_dict[test][1], True)
        result = interval1.intersection(interval2).as_list_with_included_ub()
        print(f"result : {result}\n")
        print(test_dict[test][2])
        assert result == test_dict[test][2]
        
        
# test function for the 'intervals_utils.py' class method 'union'
def test_union():
    
    # dictionary of inputs and expected ouputs for the method
    test_dict = {
        
        "identical" : [[100, 129, 150, 209, 240, 299],
                       [100, 129, 150, 209, 240, 299],
                   [100, 129, 150, 209, 240, 299]],
        
        "minus-CDS" : [[100, 129, 150, 209, 240, 299],
                       [100, 129, 240, 299],
                   [100, 129, 150, 209, 240, 299]],
        
        "fusion" : [[100, 129, 150, 209, 240, 299],
                    [100, 209, 240, 299],
                   [100, 209, 240, 299]],
        
        "shift" : [[100, 129, 150, 209, 240, 299],
                   [100, 129, 151, 209, 240, 299],
                   [100, 129, 150, 209, 240, 299]],
        
        "reverse" : [[0, 59, 90, 149, 170, 199],
                     [0, 59, 90, 149, 170, 199],
                     [0, 59, 90, 149, 170, 199]],
        
        "reverse-modified" : [[0, 59, 90, 149, 170, 199],
                              [0, 59, 80, 149, 170, 199],
                              [0, 59, 80, 149, 170, 199]],
        
        "diff-start-before" : [[100, 129, 150, 209, 240, 299],
                               [40, 69, 90, 149, 180, 239],
                               [40, 69, 90, 299]],
        
        "diff-start-after" : [[100, 129, 150, 209, 240, 299],
                              [160, 189, 210, 269, 300, 359],
                              [100, 129, 150, 359]],
        
        }
    
    print("\n*************Testing the intersection function*************")
    
    for test in test_dict:
        
        print(f"\n{test} file test")
        
        interval1 = iu.OrderedIntervals(test_dict[test][0], True)
        interval2 = iu.OrderedIntervals(test_dict[test][1], True)
        result = interval1.union(interval2).as_list_with_included_ub()
        print(f"result : {result}\n")
        print(test_dict[test][2])
        assert result == test_dict[test][2]
        

