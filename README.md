# Bézier Interpolation

# Problem
- Input: points
- Output: control points of a Bézier that go through the input points

# Algorithm
## Method 1(all in matrix):
1. Determine the initial set of nodes.
2. Solve the linear least squares parametric functional problem for the control points.
3. Solve the nonlinear least squares problem for the nearest points(this improve the curve using Newton's method).
4. Repeat steps two and three until the algorithm reaches the stopping criteria.
## Medthod 2: reuse of method 1, but instead of use 1 curve, this method use multi degree 3 Bézier curve

# Result

## Test case 1(file 1.txt):
- Grey line is the first result, black line is the last result after the line has improved. 
![Test case 1 with degree 4 Bézier curve](https://imgur.com/crgT8Y8)
![Test case 1 with multi degree 3 Bézier curve](https://imgur.com/fOQ2zwh)
## Test case 2(file 2.txt):
![Test case 2 with degree 8 Bézier curve](https://imgur.com/q90lH2o)
![Test case 2 with multi degree 3 Bézier curve](https://imgur.com/oSrN1KP)
## Test case 3(file 3.txt):
![Test case 3 with degree 5 Bézier curve](https://imgur.com/pNGi9U7)
![Test case 3 with multi degree 3 Bézier curve](https://imgur.com/5RozUNL)
## Test case 4(file 4.txt):
![Test case 4 with degree 7 Bézier curve](https://imgur.com/vLq6QzE)
![Test case 4 with multi degree 3 Bézier curve](https://imgur.com/GecM3yt)