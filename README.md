# XY-sorting
X + Y sorting in O(n ^ 2) with CDF approximations.

There are a lot of unsolved problems out there, under which: https://en.wikipedia.org/wiki/X_%2B_Y_sorting.
By merging different ideas (*) into a robust approach, we attain a good-enough approximation for the CDF of the sorted X + Y.

(*) - The ideas are mainly from:
1) https://arxiv.org/abs/1712.01208 - "to speed-up sorting is to use an existing CDF model to put the records roughly in sorted order and then correct the nearly perfectly sorted data, for example, with insertion sort"
2) http://databasearchitects.blogspot.com/2019/05/why-use-learning-when-you-can-fit.html - using a spline to fit the CDF.

Milestones:
1) Try different variants of the approximation of CDF by using only n spline knots, i.e, sqrt(size(X + Y)).
2) Decide to keep only the artificial method - take all the sums of type X[index] + Y[index] as support points for CDF.
3) Try to supress the duplicates in the CDF of X + Y, by using a hash-table. In this way, we obtain the final order faster and in a more comfortable way. 
