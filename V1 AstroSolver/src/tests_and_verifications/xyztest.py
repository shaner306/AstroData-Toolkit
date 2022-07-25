def merge_sorted_array(nums1,nums2):
    i = 0
    j = 0
    while i < len(nums1) and j < len(nums2):
        if nums1[i] < nums2[j]:
            i += 1
        else:
            nums1.insert(i,nums2[j])
            j += 1
    if j < len(nums2):
        nums1[i:] = nums2[j:]
    return nums1