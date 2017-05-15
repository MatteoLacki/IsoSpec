/*
 * MIT License
 *
 * Copyright (c) 2016 Luca Weihs
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * Implements the RangeTree datastructure.
 *
 * See the documentation of the class RangeTree for more details regarding the purpose
 * of RangeTree's.
 *
 * This file contains two main abstractions meant for use:
 * 1. The RangeTree class.
 * 2. A Point class which captures the idea of a d-dimensional euclidean point.
 */

#ifndef RANGETREE_H
#define RANGETREE_H

#include <vector>
#include <iostream>
#include <sstream>
#include <numeric>

namespace RangeTree {

    /**
     * A point in euclidean space.
     *
     * A class that represents a multi-dimensional euclidean point
     * with some associated value. We allow for each point to have an
     * associated value so that some more information can be stored with
     * each point. Points can also have a multiplicity/count, this corresponds
     * to having several duplicates of the same point.
     */
    template<class T>
    class Point {
    private:
        std::vector<double> vec;
        T val;
        int multiplicity;

    public:
        /**
         * Constructs an empty point.
         *
         * Creates a point in 0 dimensional euclidean space. This constructor
         * is provided only to make certain edge cases easier to handle.
         */
        Point() : multiplicity(0) {}

        /**
         * Constructs a point.
         *
         * Creates a point with its position in euclidean space defined by vec,
         * value defined by val, and a multiplicity/count of 1.
         *
         * @param vec the position in euclidean space.
         * @param val the value associated with the point.
         */
        Point(const std::vector<double>& vec, const T& val): val(val), vec(vec), multiplicity(1) {}

        /**
         * Constructs a point.
         *
         * Copies a point.
         *
         * @param vec the position in euclidean space.
         * @param val the value associated with the point.
         */
        Point(const Point<T>& p): val(p.val), vec(p.vec), multiplicity(p.count()) {}


        /**
         * Euclidean position of the point.
         *
         * @return the euclidean position of the point as a std::vector.
         */
        const std::vector<double>& asVector() const {
            return vec;
        }

        /**
         * The point's ambient dimension.
         *
         * @return the dimension of the space in which the point lives. I.e. a point of the
         *         form (1,2,3) lives in dimension 3.
         */
        unsigned long dim() const {
            return vec.size();
        }

        /**
         * The point's count/multiplicity.
         *
         * @return returns the count/multiplicity.
         */
        int count() const {
            return multiplicity;
        }

        /**
         * Increase the point's count/multiplicity.
         *
         * @param n amount to increase by.
         */
        void increaseCountBy(int n) {
            if (n < 0) {
                throw std::logic_error("Can't increase by a negative amount");
            }
            multiplicity += n;
        }

        /**
         * The point's value.
         *
         * @return the value stored in the point.
         */
        T value() const {
            return val;
        }

        /**
         * Index a point.
         *
         * Get the ith coordinate value of the point. I.e. if a point is of the form (4,5,6),
         * then its 0th coordinate value is 4 while its 2nd is 6.
         *
         * @param index the coordinate to index.
         * @return the coordinate value.
         */
        double operator[](int index) const {
            if(index < 0 || index >= dim()) {
                throw std::out_of_range("[] access index for point is out of range.");
            }
            return vec[index];
        }

        /**
         * Check for equality.
         *
         * Two points are considered equal if they are in the same spot, have the same
         * multiplicity/count, and store the same value.
         *
         * @param p some other point
         * @return true if \p equals the current point, otherwise false.
         */
        double operator==(const Point<T>& p) const {
            return vec == p.vec && multiplicity == p.multiplicity && val == p.val;
        }

        /**
         * Check for inequality.
         *
         * The opposite of ==.
         *
         * @param p some other point.
         * @return false if \p equals the current point, otherwise true.
         */
        double operator!=(const Point<T>& p) const {
            return !((*this) == p);
        }

        /**
         * Prints the point to standard out.
         *
         * As an example, a point with euclidean location (3,4,5) and with a
         * multiplicity/count of 4 will be printed as
         *
         * (3, 4, 5) : 4
         *
         * @param withCount whether or not to display the points count/multiplicity.
         */
        void print(bool withCount=true) const {
            std::cout << "(";
            for (int i = 0; i < dim() - 1; i++) {
                std::cout << (*this)[i] << ", ";
            }
            if (withCount) {
                std::cout << (*this)[dim() - 1] << ") : " << count() << std::endl;
            } else {
                std::cout << (*this)[dim() - 1] << ") : " << std::endl;
            }
        }
    };

    /**
     * A class that totally orders Point<T>'s in euclidean space.
     *
     * A total order of Points is required in the RangeTree. This is an implementation
     * detail that can be ignored. Given a start index \compareStartIndex, this class
     * orders points so that, for two points p_1 = (p_{11}, p_{12},...,p_{1n}) and
     * p_2 = (p_{21}, p_{22},...,p_{2n}) we have that p_1 < p_2 if and only if
     *
     * (p_{1\compareStartInd},...,p_{1n}) < (p_{2\compareStartInd},...,p_{2n})
     *
     * using the usual lexicographic order, or
     *
     * (p_{1\compareStartInd},...,p_{1n}) == (p_{2\compareStartInd},...,p_{2n}) and
     * (p_{11},...,p_{1(\compareStartInd-1)}) < (p_{21},...,p_{2(\compareStartInd-1)})
     *
     * again using the usual lexicographic order.
     */
    template <class T>
    class PointOrdering {
    private:
        int compareStartIndex;

    public:
        PointOrdering(int compareStartIndex): compareStartIndex(compareStartIndex) {
            if (compareStartIndex < 0) {
                throw new std::logic_error("Cannot have comparison start index <0.");
            }
        }

        static bool equals(const Point<T>& p1, const Point<T>& p2) {
            return p1.asVector() == p2.asVector();
        }

        int getCompareStartIndex() const {
            return compareStartIndex;
        }

        bool less(const Point<T>& p1, const Point<T>& p2) const {
            if (p1.dim() != p2.dim()) {
                throw std::logic_error("Points are incomparable (differing dims).");
            }
            if (compareStartIndex >= p1.dim()) {
                throw std::logic_error("Cannot compare points, compare start index >= point dimension.");
            }
            for (int i = compareStartIndex; i < p1.dim(); i++) {
                if (p1[i] < p2[i]) {
                    return true;
                } else if (p1[i] > p2[i]) {
                    return false;
                }
            }
            for (int i = 0; i < compareStartIndex; i++) {
                if (p1[i] < p2[i]) {
                    return true;
                } else if (p1[i] > p2[i]) {
                    return false;
                }
            }
            return false;
        }

        bool lessOrEq(const Point<T>& p1, const Point<T>& p2) const {
            return less(p1, p2) || equals(p1, p2);
        }

        bool greater(const Point<T>& p1, const Point<T>& p2) const {
            return less(p2, p1);
        }

        bool greaterOrEq(const Point<T>& p1, const Point<T>& p2) const {
            return greater(p1, p2) || equals(p1,p2);
        }

        bool operator()(const Point<T>& p1, const Point<T>& p2) const {
            return this->less(p1, p2);
        }
    };

    /**
     * A class representing a single node in a RangeTree. These should not be
     * constructed directly, instead use the RangeTree class.
     */
    template <class T>
    class RangeTreeNode {
    private:
        std::shared_ptr<RangeTreeNode<T> > left; /**< Contains points <= the comparison point **/
        std::shared_ptr<RangeTreeNode<T> > right; /**< Contains points > the comparison point **/
        std::shared_ptr<RangeTreeNode<T> > treeOnNextDim; /**< Tree on the next dimension **/
        Point<T>* point; /**< The comparison point **/
        bool isLeaf; /**< Whether or not the point is a leaf **/
        int pointCountSum; /**< Total number of points, counting multiplicities, at leaves of the tree **/
        PointOrdering<T> pointOrdering; /**< Helper to totally order input points **/

        /**
         * Create a range tree structure from input points.
         *
         * Let P = {p_1,...,p_n} be a collection of Point<T>'s that have be sorted according to the
         * lexicographic order defined by \compareStartInd (see PointOrdering) and which contain no duplicates
         * (i.e. p_i is strictly less than p_{i+1} for all i). Then creates a range tree structure on the points
         * p_{first} to p_{last} given the lexicographic order (when appropriate, sub range tree structures are
         * created using a lexicographic order defined by \compareStartInd + 1).
         *
         * @param sortedUniquePoints a std::vector of Point<T>s sorted according to the lexicographic order defined
         *                           above (these points must be unique with respect to that order).
         * @param first an int defining the first point to use from the input collection.
         * @param last an int defining the last point to use from the input collection.
         * @param compareStartInd an int defining the PointOrdering lexicographic order.
         * @return a RangeTreeNode representing the root of a new range tree structure.
         */
        std::shared_ptr<RangeTreeNode<T> > sortedPointsToBinaryTree(
                const std::vector<Point<T>* >& sortedUniquePoints,
                int first, int last, int compareStartInd) {
            if (sortedUniquePoints.size() == 0) {
                throw std::logic_error("Number of points input to sortedPointsToBinaryTree must be >0.");
            }

            if (first > last || first < 0 || last >= sortedUniquePoints.size()) {
                std::ostringstream out;
                out << "Must have 0 <= first (" << first << ")  <= last (" << last << "), and "
                        "last < the number of input points (" << sortedUniquePoints.size() << ").";
                throw std::logic_error(out.str());
            } else if (first == last) {
                return std::shared_ptr<RangeTreeNode<T>>(
                        new RangeTreeNode(sortedUniquePoints[first], compareStartInd));
            } else {
                int mid = (first + last) / 2;
                auto left = sortedPointsToBinaryTree(sortedUniquePoints, first, mid, compareStartInd);
                auto right = sortedPointsToBinaryTree(sortedUniquePoints, mid + 1, last, compareStartInd);

                std::vector<Point<T>* > subVec;
                for (int i = first; i <= last; i++) {
                    subVec.push_back(sortedUniquePoints[i]);
                }

                return std::shared_ptr<RangeTreeNode<T>>(
                        new RangeTreeNode(left, right,
                                          subVec,
                                          sortedUniquePoints[mid],
                                          compareStartInd));
            }
        }

    public:
        /**
         * Construct a range tree structure from points.
         *
         * Creates a range tree structure on the input collection \allPoints using the lexicographic order
         * starting at \compareStartInd.
         *
         * @param uniquePoints a collection of points.
         * @param compareStartInd the index to use for the lexicographic order
         * @return a range tree structure
         */
        RangeTreeNode(std::vector<Point<T>* > uniquePoints, int compareStartInd): pointOrdering(compareStartInd) {
            if (uniquePoints.size() == 0) {
                throw std::range_error("Range tree requires input vector of points to not be empty.");
            }

//            sort(idx.begin(), idx.end(),
//             [&points, &pointOrdering](size_t i1, size_t i2) {
//                 return pointOrdering.less(points[i1], points[i2]); });

            std::sort(uniquePoints.begin(), uniquePoints.end(),
                      [this](const Point<T>* p1, const Point<T>* p2) {
                return this->pointOrdering.less(*p1, *p2);
            });

            int numLeafNodes = uniquePoints.size();
            int mid = (numLeafNodes - 1) / 2;
            point = uniquePoints[mid];

            if (numLeafNodes == 1) {
                isLeaf = true;
                pointCountSum = point->count();
            } else {
                left = sortedPointsToBinaryTree(uniquePoints, 0, mid, compareStartInd);
                right = sortedPointsToBinaryTree(uniquePoints, mid + 1, numLeafNodes - 1, compareStartInd);
                pointCountSum = left->totalPoints() + right->totalPoints();
                if (compareStartInd + 1 != point->dim()) {
                    treeOnNextDim = std::shared_ptr<RangeTreeNode>(
                            new RangeTreeNode(uniquePoints, compareStartInd + 1));
                }
                isLeaf = false;
            }
        }

        /**
         * Construct a range tree structure
         *
         * Creates a range tree structure similarily as for
         * RangeTreeNode(std::vector<Point<T> > sortedUniquePoints, int compareStartInd) but where the root node's left and
         * right subtrees are already known.
         *
         * @param left std::shared_ptr to left subtree
         * @param right std::shared_ptr to right subtree
         * @param sortedUniquePoints collection of sorted unqiue points
         * @param comparePoint the point at the node to be constructed that is used as the comparison point
         * @param compareStartInd the index defining the lexicographic order
         */
        RangeTreeNode(const std::shared_ptr<RangeTreeNode<T> >& left,
                      const std::shared_ptr<RangeTreeNode<T> >& right,
                      const std::vector<Point<T>* >& sortedUniquePoints,
                      Point<T>* comparePoint,
                      int compareStartInd) :
                left(left), right(right),
                point(comparePoint), isLeaf(false),
                pointOrdering(compareStartInd) {
            if (compareStartInd + 1 != point->dim()) {
                treeOnNextDim = std::shared_ptr<RangeTreeNode<T>>(new RangeTreeNode<T>(sortedUniquePoints,
                                                                                       compareStartInd + 1));
            }
            pointCountSum = (*left).totalPoints() + (*right).totalPoints();
        }

        /**
         * Construct a RangeTreeNode representing a leaf.
         *
         * @param pointAtLeaf the point to use at the leaf.
         * @param compareStartInd the index defining the lexicographic ordering.
         * @return
         */
        RangeTreeNode(Point<T>* pointAtLeaf, int compareStartInd) :
                point(pointAtLeaf), isLeaf(true), pointCountSum(pointAtLeaf->count()), pointOrdering(compareStartInd) {}

        /**
         * Total count of points at the leaves of the range tree rooted at this node.
         *
         * The total count returned INCLUDES the multiplicities/count of the points at the leaves.
         *
         * @return the total count.
         */
        int totalPoints() const {
            return pointCountSum;
        }

        /**
         * Return all points at the leaves of the range tree rooted at this node.
         * @return all the points.
         */
        std::vector<Point<T> > getAllPoints() const {
            if (isLeaf) {
                std::vector<Point<T> > vec;
                vec.push_back(*point);
                return vec;
            }
            auto allPointsLeft = (*left).getAllPoints();
            auto allPointsRight = (*right).getAllPoints();

            allPointsLeft.insert(allPointsLeft.end(), allPointsRight.begin(), allPointsRight.end());
            return allPointsLeft;
        }

        /**
         * Check if point is in a euclidean box.
         *
         * Determines whether or not a point is in a euclidean box. That is, if p_1 = (p_{11},...,p_{1n}) is an
         * n-dimensional point. Then this function returns true if, for all 1 <= i <= n we have
         *
         * lower[i] <= p_{1i} <= upper[i] if withLower[i] == true and withUpper[i] == true, or
         * lower[i] < p_{1i} <= upper[i] if withLower[i] == false and withUpper[i] == true, or
         * lower[i] <= p_{1i} < upper[i] if withLower[i] == true and withUpper[i] == false, or
         * lower[i] < p_{1i} < upper[i] if withLower[i] == false and withUpper[i] == false.
         *
         * @param point the point to check.
         * @param lower the lower points of the rectangle.
         * @param upper the upper bounds of the rectangle.
         * @param withLower whether to use strict (<) or not strict (<=) inequalities at certain coordiantes of p_1
         *                  for the lower bounds.
         * @param withUpper as for \withLower but for the upper bounds.
         * @return true if the point is in the rectangle, false otherwise.
         */
        bool pointInRange(const Point<T>& point,
                          const std::vector<double>& lower,
                          const std::vector<double>& upper,
                          const std::vector<bool>& withLower,
                          const std::vector<bool>& withUpper) const {
            for (int i = 0; i < point.dim(); i++) {
                if (point[i] < lower[i] ||
                    (!withLower[i] && point[i] == lower[i])) {
                    return false;
                }
                if (point[i] > upper[i] ||
                    (!withUpper[i] && point[i] == upper[i])) {
                    return false;
                }
            }
            return true;
        }

        bool pointInRange(const Point<T>* point,
                          const std::vector<double>& lower,
                          const std::vector<double>& upper,
                          const std::vector<bool>& withLower,
                          const std::vector<bool>& withUpper) const {
            return pointInRange(*point, lower, upper, withLower, withUpper);
        }

        /**
         * Count the number of points at leaves of tree rooted at the current node that are within the given bounds.
         *
         * @param lower see the pointInRange(...) function.
         * @param upper
         * @param withLower
         * @param withUpper
         * @return the count.
         */
        int countInRange(const std::vector<double>& lower,
                         const std::vector<double>& upper,
                         const std::vector<bool>& withLower,
                         const std::vector<bool>& withUpper) const {
            if (isLeaf) {
                if (pointInRange(point, lower, upper, withLower, withUpper)) {
                    return totalPoints();
                } else {
                    return 0;
                }
            }
            int compareInd = pointOrdering.getCompareStartIndex();

            if (((*point)[compareInd] > upper[compareInd]) ||
                ((*point)[compareInd] == upper[compareInd] && !withUpper[compareInd])) {
                return (*left).countInRange(lower, upper, withLower, withUpper);
            }
            if (((*point)[compareInd] < lower[compareInd]) ||
                ((*point)[compareInd] == lower[compareInd] && !withLower[compareInd])) {
                return (*right).countInRange(lower, upper, withLower, withUpper);
            }

            std::vector<std::shared_ptr<RangeTreeNode<T> > > canonicalNodes;

            if ((*left).isLeaf) {
                canonicalNodes.push_back(left);
            } else {
                (*left).leftCanonicalNodes(lower, withLower, canonicalNodes);
            }

            if ((*right).isLeaf) {
                canonicalNodes.push_back(right);
            } else {
                (*right).rightCanonicalNodes(upper, withUpper, canonicalNodes);
            }

            int numPointsInRange = 0;
            for (int i = 0; i < canonicalNodes.size(); i++) {
                std::shared_ptr<RangeTreeNode<T> > node = canonicalNodes[i];
                if ((*node).isLeaf) {
                    if (pointInRange((*node).point, lower, upper, withLower, withUpper)) {
                        numPointsInRange += (*node).totalPoints();
                    }
                } else if (compareInd + 1 == point->dim()) {
                    numPointsInRange += (*node).totalPoints();
                } else {
                    numPointsInRange += node->treeOnNextDim->countInRange(lower, upper, withLower, withUpper);
                }
            }

            return numPointsInRange;
        }

        /**
         * Return the points at leaves of tree rooted at the current node that are within the given bounds.
         *
         * @param lower see the pointInRange(...) function.
         * @param upper
         * @param withLower
         * @param withUpper
         * @return a std::vector of the Points.
         */
        std::vector<Point<T> > pointsInRange(const std::vector<double>& lower,
                                             const std::vector<double>& upper,
                                             const std::vector<bool>& withLower,
                                             const std::vector<bool>& withUpper) const {
            std::vector<Point<T> > pointsToReturn = {};
            if (isLeaf) {
                if (pointInRange(point, lower, upper, withLower, withUpper)) {
                    pointsToReturn.push_back(*point);
                }
                return pointsToReturn;
            }
            int compareInd = pointOrdering.getCompareStartIndex();

            if (((*point)[compareInd] > upper[compareInd]) ||
                ((*point)[compareInd] == upper[compareInd] && !withUpper[compareInd])) {
                return (*left).pointsInRange(lower, upper, withLower, withUpper);
            }
            if (((*point)[compareInd] < lower[compareInd]) ||
                ((*point)[compareInd] == lower[compareInd] && !withLower[compareInd])) {
                return (*right).pointsInRange(lower, upper, withLower, withUpper);
            }

            std::vector<std::shared_ptr<RangeTreeNode<T> > > canonicalNodes = {};

            if ((*left).isLeaf) {
                canonicalNodes.push_back(left);
            } else {
                (*left).leftCanonicalNodes(lower, withLower, canonicalNodes);
            }

            if ((*right).isLeaf) {
                canonicalNodes.push_back(right);
            } else {
                (*right).rightCanonicalNodes(upper, withUpper, canonicalNodes);
            }

            for (int i = 0; i < canonicalNodes.size(); i++) {
                std::shared_ptr<RangeTreeNode<T> > node = canonicalNodes[i];
                if ((*node).isLeaf) {
                    if (pointInRange((*node).point, lower, upper, withLower, withUpper)) {
                        pointsToReturn.push_back(*(node->point));
                    }
                } else if (compareInd + 1 == point->dim()) {
                    auto allPointsAtNode = node->getAllPoints();
                    pointsToReturn.insert(pointsToReturn.end(), allPointsAtNode.begin(), allPointsAtNode.end());
                } else {
                    auto allPointsAtNode = node->treeOnNextDim->pointsInRange(lower, upper, withLower, withUpper);
                    pointsToReturn.insert(pointsToReturn.end(), allPointsAtNode.begin(), allPointsAtNode.end());
                }
            }

            return pointsToReturn;
        }

        /**
         * Helper function for countInRange(...).
         * @param lower
         * @param withLower
         * @param nodes
         */
        void leftCanonicalNodes(const std::vector<double>& lower,
                                const std::vector<bool>& withLower,
                                std::vector<std::shared_ptr<RangeTreeNode<T> > >& nodes) {
            if (isLeaf) {
                throw std::logic_error("Should never have a leaf deciding if its canonical.");
            }
            int compareInd = pointOrdering.getCompareStartIndex();
            int totalPoints = 0;
            if (lower[compareInd] < (*point)[compareInd] ||
                (lower[compareInd] == (*point)[compareInd] && withLower[compareInd])) {
                nodes.push_back(right);
                if ((*left).isLeaf) {
                    nodes.push_back(left);
                } else {
                    (*left).leftCanonicalNodes(lower, withLower, nodes);
                }
            } else {
                if ((*right).isLeaf) {
                    nodes.push_back(right);
                } else {
                    (*right).leftCanonicalNodes(lower, withLower, nodes);
                }
            }
        }

        /**
         * Helper function for countInRange(...).
         * @param upper
         * @param withUpper
         * @param nodes
         */
        void rightCanonicalNodes(const std::vector<double>& upper,
                                 const std::vector<bool>& withUpper,
                                 std::vector<std::shared_ptr<RangeTreeNode<T> > >& nodes) {
            if (isLeaf) {
                throw std::logic_error("Should never have a leaf deciding if its canonical.");
            }
            int compareInd = pointOrdering.getCompareStartIndex();
            int totalPoints = 0;
            if (upper[compareInd] > (*point)[compareInd] ||
                (upper[compareInd] == (*point)[compareInd] && withUpper[compareInd])) {
                nodes.push_back(left);
                if ((*right).isLeaf) {
                    nodes.push_back(right);
                } else {
                    (*right).rightCanonicalNodes(upper, withUpper, nodes);
                }
            } else {
                if ((*left).isLeaf) {
                    nodes.push_back(left);
                } else {
                    (*left).rightCanonicalNodes(upper, withUpper, nodes);
                }
            }
        }

        /**
         * Print the structure of the tree rooted at the curret node.
         *
         * The printed structure does not reflect any subtrees for other coordinates.
         *
         * @param numIndents the number of indents to use before every line printed.
         */
        void print(int numIndents) {
            for (int i = 0; i < numIndents; i++) { std::cout << "\t"; }
            if (isLeaf) {
                point->print(true);
            } else {
                point->print(false);
                left->print(numIndents + 1);
                right->print(numIndents + 1);
            }
        }
    };

    /**
     * A class facilitating fast orthogonal range queries.
     *
     * A RangeTree allows for 'orthogonal range queries.' That is, given a collection of
     * points P = {p_1, ..., p_n} in euclidean d-dimensional space, a RangeTree can efficiently
     * answer questions of the form
     *
     * "How many points of p are in the box high dimensional rectangle
     * [l_1, u_1] x [l_2, u_2] x ... x [l_d, u_d]
     * where l_1 <= u_1, ..., l_n <= u_n?"
     *
     * It returns the number of such points in worst case
     * O(log(n)^d) time. It can also return the points that are in the rectangle in worst case
     * O(log(n)^d + k) time where k is the number of points that lie in the rectangle.
     *
     * The particular algorithm implemented here is described in Chapter 5 of the book
     *
     * Mark de Berg, Otfried Cheong, Marc van Kreveld, and Mark Overmars. 2008.
     * Computational Geometry: Algorithms and Applications (3rd ed. ed.). TELOS, Santa Clara, CA, USA.
     */
    template <class T>
    class RangeTree {
    private:
        std::shared_ptr<RangeTreeNode<T>> root;
        std::vector<std::shared_ptr<Point<T> > > sortedUniquePoints;

    public:
        /**
         * Construct a new RangeTree from input points.
         *
         * Input points may have duplicates but if two points p_1 and p_2 are at the same euclidean position
         * then they are required to have the same value as all duplicate points will be accumulated into a
         * single point with multiplicity/count equal to the sum of the multiplicities/counts of all such duplicates.
         *
         * @param points the points from which to create a RangeTree
         */
        RangeTree(const std::vector<Point<T> >& points) {
            if (points.size() != 0) {
                int dim = points[0].dim();
                for (int i = 1; i < points.size(); i++) {
                    if (points[i].dim() != dim) {
                        throw std::logic_error("Input points to RangeTree must all have the same dimension.");
                    }
                }
            }

            PointOrdering<T> pointOrdering(0);

            sortedUniquePoints.push_back(std::shared_ptr<Point<T> >(new Point<T>(points[0])));
            int k = 0;
            for (int i = 1; i < points.size(); i++) {
                if (pointOrdering.equals(*(sortedUniquePoints[k]), points[i])) {
                    if (sortedUniquePoints[k]->value() != points[i].value()) {
                        throw std::logic_error("Input points have same position but different values");
                    }
                    sortedUniquePoints[k]->increaseCountBy(points[i].count());
                } else {
                    sortedUniquePoints.push_back(std::shared_ptr<Point<T> >(new Point<T>(points[i])));
                    k++;
                }
            }

            std::vector<Point<T>* > sortedUniquePointsRaw;
            for (int i = 0; i < sortedUniquePoints.size(); i++) {
                sortedUniquePointsRaw.push_back(&(*sortedUniquePoints[i]));
            }

            root = std::shared_ptr<RangeTreeNode<T>>(new RangeTreeNode<T>(sortedUniquePointsRaw, 0));
        }

        /**
         * The number of points within a high dimensional rectangle.
         *
         * The rectangle is defined by the input parameters. In particular, an n-dimensional point
         * p_1 = (p_{11},...,p_{1n}) is an is in the rectangle if, for all 1 <= i <= n, we have
         *
         * lower[i] <= p_{1i} <= upper[i] if withLower[i] == true and withUpper[i] == true, or
         * lower[i] < p_{1i} <= upper[i] if withLower[i] == false and withUpper[i] == true, or
         * lower[i] <= p_{1i} < upper[i] if withLower[i] == true and withUpper[i] == false, or
         * lower[i] < p_{1i} < upper[i] if withLower[i] == false and withUpper[i] == false.
         *
         * @param lower the lower bounds of the rectangle.
         * @param upper the upper bounds of the rectangle.
         * @param withLower whether to use strict (<) or not strict (<=) inequalities at certain coordiantes of p_1
         *                  for the lower bounds.
         * @param withUpper as for \withLower but for the upper bounds.
         * @return the number of points in the rectangle.
         */
        int countInRange(const std::vector<double>& lower,
                         const std::vector<double>& upper,
                         const std::vector<bool>& withLower,
                         const std::vector<bool>& withUpper) const {
            return (*root).countInRange(lower, upper, withLower, withUpper);
        }

        /**
         * Return all points in range.
         *
         * Returns a std::vector of all points in the given rectangle. See \countInRange for how
         * this rectangle is specified. NOTE: these points may not be identical to those points
         * that were given as input to the RangeTree at construction time. This is because
         * duplicate points are merged together with appropriate incrementing of their multiplicity.
         * That is, two points at euclidean position (1,2,3) and multiplicities/counts of 2 and 3
         * respectively will be merged into a single Point with position (1,2,3) and multiplicity 5
         * (recall that all points with the same euclidean position are required to have the same
         * associated value so it is ok to merge in this way).
         *
         * @param lower the lower bounds of the rectangle.
         * @param upper the upper bounds of the rectangle.
         * @param withLower whether to use strict (<) or not strict (<=) inequalities at certain coordiantes of p_1
         *                  for the lower bounds.
         * @param withUpper as for \withLower but for the upper bounds.
         * @return the number of points in the rectangle.
         */
        std::vector<Point<T> > pointsInRange(const std::vector<double>& lower,
                                             const std::vector<double>& upper,
                                             const std::vector<bool>& withLower,
                                             const std::vector<bool>& withUpper) const {
            return (*root).pointsInRange(lower, upper, withLower, withUpper);
        }

        void print() const {
            root->print(0);
        }
    };

    /**
     * A class which is used to naively count the number of points in a given rectangle. This class is used
     * for testing an benchmarking, it should not be used in practice.
     */
    template <class T>
    class NaiveRangeCounter {
    private:
        std::vector<Point<T> > points;

        static bool pointInRange(const Point<T>& point,
                                 const std::vector<double>& lower,
                                 const std::vector<double>& upper,
                                 const std::vector<bool>& withLower,
                                 const std::vector<bool>& withUpper) {
            for (int i = 0; i < point.dim(); i++) {
                if (point[i] < lower[i] ||
                    (point[i] == lower[i] && !withLower[i])) {
                    return false;
                }
                if (point[i] > upper[i] ||
                    (point[i] == upper[i] && !withUpper[i])) {
                    return false;
                }
            }
            return true;
        }

    public:
        NaiveRangeCounter(std::vector<Point<T> > points): points(points) {}

        int countInRange(const std::vector<double>& lower,
                         const std::vector<double>& upper,
                         const std::vector<bool>& withLower,
                         const std::vector<bool>& withUpper) const {
            int count = 0;
            for (int i = 0; i < points.size(); i++) {
                if (pointInRange(points[i], lower, upper, withLower, withUpper)) {
                    count += points[i].count();
                }
            }
            return count;
        }

        std::vector<Point<T> > pointsInRange(const std::vector<double>& lower,
                                             const std::vector<double>& upper,
                                             const std::vector<bool>& withLower,
                                             const std::vector<bool>& withUpper) const {
            std::vector<Point<T> > selectedPoints = {};
            for (int i = 0; i < points.size(); i++) {
                if (pointInRange(points[i], lower, upper, withLower, withUpper)) {
                    selectedPoints.push_back(points[i]);
                }
            }
            return selectedPoints;
        }
    };

} // namespace

#endif //RANGETREE_H
