# Bézier Interpolation

## Problem
- Input: points
- Output: control points of a Bézier that go through the input points

## Algorithm
### Method 1(all in matrix):
1. Determine the initial set of nodes.
2. Solve the linear least squares parametric functional problem for the control points.
3. Solve the nonlinear least squares problem for the nearest points(this improve the curve using Newton's method).
4. Repeat steps two and three until the algorithm reaches the stopping criteria.
### Medthod 2: reuse of method 1, but instead of use 1 curve, this method use multi degree 3 Bézier curve

## Result

### P/S
- This program allows 2 and 3 dimension, that is, case 1 and 2 use 3D points, and to make it easier to observe, the result is looked from Oxy and Oyz plan.
- Grey line is the first result, black line is the last result after the line has improved. 
- In each test case, the first image is the result of the method 1, the second image is the result of the method 2.
### Test case 1(file 1.txt):
- Test case 1 with degree 4 Bézier curve
![Test case 1 with degree 4 Bézier curve](https://i.imgur.com/crgT8Y8.png)
- Test case 1 with multi degree 3 Bézier curve
![Test case 1 with multi degree 3 Bézier curve](https://i.imgur.com/fOQ2zwh.png)
### Test case 2(file 2.txt):
- Test case 2 with degree 8 Bézier curve
![Test case 2 with degree 8 Bézier curve](https://i.imgur.com/q90lH2o.png)
- Test case 2 with multi degree 3 Bézier curve
![Test case 2 with multi degree 3 Bézier curve](https://i.imgur.com/oSrN1KP.png)
### Test case 3(file 3.txt):
- Test case 3 with degree 5 Bézier curve
![Test case 3 with degree 5 Bézier curve](https://i.imgur.com/pNGi9U7.png)
- Test case 3 with multi degree 3 Bézier curve
![Test case 3 with multi degree 3 Bézier curve](https://i.imgur.com/5RozUNL.png)
### Test case 4(file 4.txt):
- Test case 4 with degree 7 Bézier curve
![Test case 4 with degree 7 Bézier curve](https://i.imgur.com/vLq6QzE.png)
- Test case 4 with multi degree 3 Bézier curve
![Test case 4 with multi degree 3 Bézier curve](https://i.imgur.com/GecM3yt.png)