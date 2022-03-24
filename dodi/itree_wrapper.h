#include <cstdint>
#include <iostream>
#include <vector>


#include "IITreeBFS.h"


struct Interval { int low, high; };

class BasicIntervalTree
// https://www.geeksforgeeks.org/interval-tree/
{
    public:

        BasicIntervalTree() {}
        ~BasicIntervalTree() {}

        void add(int start, int end, int index) {
            tree.add(start, end, index);
        }

        void index() {
            tree.index();
        }

        bool searchInterval(int pos, int pos2) {
            std::vector<size_t> a;
            tree.overlap(pos, pos2, a);
            if (a.size() == 0)
                return false;
            else
                return true;
        }

        void allOverlappingIntervals(int pos, int pos2, std::vector<int>& res) {
            std::vector<size_t> a;
            tree.overlap(pos, pos2, a);
            for (size_t i = 0; i < a.size(); ++i) {
                res.push_back(tree.data(a[i]));
            };
        }

        int countOverlappingIntervals(int pos, int pos2) {
            std::vector<size_t> a;
            tree.overlap(pos, pos2, a);
            return (int)a.size();
        }

        Interval* overlappingInterval(int pos, int pos2) {
            std::vector<size_t> a;
            tree.overlap(pos, pos2, a);
            Interval *res = new Interval;
            for (size_t i = 0; i < a.size(); ++i) {
                res->low = tree.start(a[i]);
                res->high = tree.end(a[i]);
                break;
            };
            return res;
        }

    private:
        IITree<int, int> tree;

};