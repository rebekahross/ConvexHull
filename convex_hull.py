from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF, QObject
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF, QObject
elif PYQT_VER == 'PYQT6':
    from PyQt6.QtCore import QLineF, QPointF, QObject
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import math
from math import atan2

# Some global color constants that might be useful
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

# Global variable that controls the speed of the recursion automation, in seconds
PAUSE = 0.25


#
# This is the class you have to complete.
#

class ConvexHullSolver:
    @staticmethod
    def divide_and_conquer(ordered_points):                                     # O(n log n) - from adding results below
        length_points = len(ordered_points)                                     # O(1) - arithmetic with fixed-length numbers
        if length_points <= 2:                                                  # O(1) - comparisons are constant time
            return ordered_points
        mid_point = length_points // 2                                          # O(1) - arithmetic with fixed-length numbers
        L = ConvexHullSolver.divide_and_conquer(ordered_points[:mid_point])     # O(log n) - n = # points
        R = ConvexHullSolver.divide_and_conquer(ordered_points[mid_point:])     # ^ Recursively calls function so half the number of points each time
        return ConvexHullSolver.merge(L, R)                                     # O(n) - see below for specifics of merge function

    @staticmethod
    def merge(L, R):             # O(2n) = O(n) (for n = # points). see while loops for more details
        length_L = len(L)        # O(1) - arithmetic with fixed-length numbers
        length_R = len(R)        # O(1) - arithmetic with fixed-length numbers
        up_left, up_right = ConvexHullSolver.upperTangent(L, length_L, R, length_R)         # O(n) - n = # points
        low_left, low_right = ConvexHullSolver.lowerTangent(L, length_L, R, length_R)       # O(n) - see below for specifics of upper/lowerTangent functions
        convex_hull = []
        i = low_left
        while True:                         # Worst case: O(n) (for n = # points), Best case: O(1)
            convex_hull.append(L[i])        # O(1) - appending to a list is constant time
            if i == up_left:                # O(1) - comparisons are constant time
                break
            i = (i + 1) % length_L          # O(1) - arithmetic of fixed-point numbers
        i = up_right
        while True:                         # Worst case: O(n) (for n = # points), Best case: O(1)
            convex_hull.append(R[i])        # O(1) - appending to a list is constant time
            if i == low_right:              # O(1) - comparisons are constant time
                break
            i = (i + 1) % length_R          # O(1) - arithmetic of fixed-point numbers
        return convex_hull

    @staticmethod
    def upperTangent(L, length_L, R, length_R):                 # O(l + r) for (l = length_L, r = length_R)
        left_index = max(range(length_L), key=lambda i: L[i].x())           # O(l) (for l = # L) because it iterates through each of the points to find the max
        right_index = min(range(length_R), key=lambda i: R[i].x())          # O(r) (for r = # R) because it iterates through each of the points to find the max
                                                                            # Best case: O(1) - finds tangent for each loop on first time through
        while True:                                                         # Worst case: O(l + r) (for l = # L and r = # R)
            done = False
            while ConvexHullSolver.matrixMultiply(R[right_index], L[left_index], L[(left_index - 1) % length_L]) > 0:       # Worst case: O(l) - has to go through all elements in L to find tangent
                left_index = (left_index - 1) % length_L                                                                    # Best case: O(1) - finds tangent on first time through
                done = True
            while ConvexHullSolver.matrixMultiply(L[left_index], R[right_index], R[(right_index + 1) % length_R]) < 0:      # Worst case: O(r) - has to go through all elements in R to find tangent
                right_index = (right_index + 1) % length_R                                                                  # Best case: O(1) - finds tangent on first time through
                done = True
            if not done:
                break
        return left_index, right_index

    @staticmethod
    def lowerTangent(L, length_L, R, length_R):                             # O(l + r) for (l = length_L, r = length_R)
        left_index = max(range(length_L), key=lambda i: L[i].x())           # O(l) (for l = # L) because it iterates through each of the points to find the max
        right_index = min(range(length_R), key=lambda i: R[i].x())          # O(r) (for r = # R) because it iterates through each of the points to find the max
                                                                            # Best case: O(1) - finds tangent for each loop on first time through
        while True:                                                         # Worst case: O(l + r) (for l = # L and r = # R)
            done = False
            while ConvexHullSolver.matrixMultiply(R[right_index], L[left_index], L[(left_index + 1) % length_L]) < 0:     # Worst case: O(l) - has to go through all elements in L to find tangent
                left_index = (left_index + 1) % length_L                                                                  # Best case: O(1) - finds tangent on first time through
                done = True
            while ConvexHullSolver.matrixMultiply(L[left_index], R[right_index], R[(right_index - 1) % length_R]) > 0:    # Worst case: O(r) - has to go through all elements in R to find tangent
                right_index = (right_index - 1) % length_R                                                                # Best case: O(1) - finds tangent on first time through
                done = True
            if not done:                        # O(1) - comparison is constant time
                break
        return left_index, right_index

    @staticmethod
    def matrixMultiply(o, a, b):                # O(1) - arithmetic function that is unaffected by the size of n
        axo = (a.x() - o.x())                   # O(1) - arithmetic of fixed-size numbers
        ayo = (a.y() - o.y())                   # O(1) - arithmetic of fixed-size numbers
        bxo = (b.x() - o.x())                   # O(1) - arithmetic of fixed-size numbers
        byo = (b.y() - o.y())                   # O(1) - arithmetic of fixed-size numbers
        return (axo * byo) - (ayo * bxo)        # O(1) - arithmetic of fixed-size numbers


    # Class constructor
    def __init__(self):
        super().__init__()
        self.pause = False

    # Some helper methods that make calls to the GUI, allowing us to send updates
    # to be displayed.

    def showTangent(self, line, color):
        self.view.addLines(line, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseTangent(self, line):
        self.view.clearLines(line)

    def blinkTangent(self, line, color):
        self.showTangent(line, color)
        self.eraseTangent(line)

    def showHull(self, polygon, color):
        self.view.addLines(polygon, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseHull(self, polygon):
        self.view.clearLines(polygon)

    def showText(self, text):
        self.view.displayStatusText(text)

    # This is the method that gets called by the GUI and actually executes
    # the finding of the hull
    def compute_hull(self, points, pause, view):
        self.pause = pause
        self.view = view
        assert (type(points) == list and type(points[0]) == QPointF)

        t1 = time.time()
        points = sorted(points, key=lambda point: point.x())
        length_points = len(points)

        t2 = time.time()

        t3 = time.time()
        points = ConvexHullSolver.divide_and_conquer(points)
        polygon = [QLineF(points[i], points[(i + 1) % (len(points))]) for i in range(len(points))]

        t4 = time.time()

        # when passing lines to the display, pass a list of QLineF objects.  Each QLineF
        # object can be created with two QPointF objects corresponding to the endpoints
        self.showHull(polygon, RED)
        self.showText('Time Elapsed (Convex Hull): {:3.3f} sec'.format(t4 - t3))

