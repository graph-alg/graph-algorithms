// The @Treap class implements the treap structure.
#pragma once
#include <vector>
#include <cctype>
#include <cstdint>

namespace gadget {
    class TreapNode {
    public:
        long long p;  // the parent
        long long l;  // the left child
        long long r;  // the right child
        long long s;  // the size of the subtree rooted at this node
        long long w;  // the priority
    };

    class Treap {
    public:
        explicit Treap(long long n);

        void Insert(long long x, bool f, long long &r);

        void InsertAfter(long long x, long long y, long long &r);

        void Delete(long long x, long long &r);

        long long Merge(long long r1, long long r2);

        long long Rank(long long x);

        long long Select(long long r, long long rank);

        long long Root(long long x);

        long long Minimum(long long x);

        long long Maximum(long long x);

        long long Size(long long r);

        void Check(long long r);

    private:

        void LeftRotate(long long &x);

        void RightRotate(long long &x);

        void SubCheck(long long x);

        long long n_;
        std::vector<TreapNode> nd_;
    };
} // namespace gadget

