/**
 * @file ExtendedUF.h
 * @brief An Union Find able to keep some data
 * @author Gueunet Charles
 * @version 1
 * @date 2016-03-30
 */

#ifndef EXTENDEDUF_H
#define EXTENDEDUF_H

#include <DataTypes.h>
#include <vector>

namespace wtfit
{
   class ExtendedUnionFind
   {
     private:
      int rank_;
      ExtendedUnionFind *parent_;
      ufDataType data_;
      idVertex origin_;

     public:
      inline ExtendedUnionFind(const idVertex &origin)
      {
         rank_      = 0;
         parent_    = this;
         data_      = nullUfData;
         origin_    = origin;
      }

      inline ExtendedUnionFind(const ExtendedUnionFind &other)
      {
         rank_      = other.rank_;
         parent_    = this;
         data_      = other.data_;
         origin_    = other.origin_;
      }

      inline void setData(const ufDataType &d)
      {
         data_ = d;
      }

      inline void setOrigin(const idVertex &origin)
      {
         origin_ = origin;
      }

      inline const ufDataType &getData(void) const
      {
         return data_;
      }

      inline const idVertex &getOrigin(void) const
      {
         return origin_;
      }

      // heavy recursif
      inline ExtendedUnionFind *find()
      {
         if (parent_ == this)
            return this;
         else {
            parent_ = parent_->find();
            return parent_;
         }
      }

      inline int getRank() const
      {
         return rank_;
      };

      inline void setParent(ExtendedUnionFind *parent)
      {
         parent_ = parent;
      };

      inline void setRank(const int &rank)
      {
         rank_ = rank;
      };

      static inline ExtendedUnionFind *makeUnion(ExtendedUnionFind *uf0, ExtendedUnionFind *uf1)
      {
         uf0 = uf0->find();
         uf1 = uf1->find();

         if (uf0 == uf1) {
            return uf0;
         } else if (uf0->getRank() > uf1->getRank()) {
            uf1->setParent(uf0);
            return uf0;
         } else if (uf0->getRank() < uf1->getRank()) {
            uf0->setParent(uf1);
            return uf1;
         } else {
            uf1->setParent(uf0);
            uf0->setRank(uf0->getRank() + 1);
            return uf0;
         }

         return NULL;
      }

      static inline ExtendedUnionFind *makeUnion(std::vector<ExtendedUnionFind *> &sets)
      {
         ExtendedUnionFind *n = NULL;

         if (!sets.size())
            return NULL;

         if (sets.size() == 1)
            return sets[0];

         for (int i = 0; i < (int)sets.size() - 1; i++)
            n = makeUnion(sets[i], sets[i + 1]);

         return n;
      }

      inline bool operator<(const ExtendedUnionFind &other) const
      {
         return rank_ < other.rank_;
      };

      inline bool operator>(const ExtendedUnionFind &other) const
      {
         return rank_ > other.rank_;
      };

   };
}

#endif /* end of include guard: EXTENDEDUF_H */

