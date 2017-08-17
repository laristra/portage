/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/



/* ---------------------------------------------------------------------------- 
 ** INCLUDES/KDTREE.H
 **
 ** Contains data structures and prototypes for manipulating polygon data.
 ** ----------------------------------------------------------------------------
 */
#ifndef _KDTREE_H_
#define _KDTREE_H_

#include <vector>
#include <set>
#include <cstdlib>
#include "BoundBox.h"
#include "portage/support/Point.h"
#include "portage/support/Vector.h"

namespace gk {

using Portage::Point;
using Portage::Vector;

#define SWAP(a, b, type) {  \
    type c_ = b; \
    b = a; \
    a = c_; \
}

const long XDIM = 0;
const long YDIM = 1;
const long ZDIM = 2;

const double KDEPS = 1.0e-8;
const double KD_SAFETY_EPS = 1.0e-4;

/*!
  @struct KDTree "kdtree.h"
  @brief An N-dimensional k-d tree for manipulating polygon data.
  @tparam D Dimension of the k-d tree.
  */
template<long D> struct KDTree {

        /// Default constructor.
        KDTree() {
            num_entities = 0;
            linkp = NULL;
            sbox = NULL;
        }
        /// Default destructor.
        ~KDTree() {
            if(linkp) {
                delete [] linkp;
            }
            if(sbox) {
                delete [] sbox;
            }
        }
        size_t num_entities;
        long *linkp;
        IsotheticBBox<D> *sbox;
};


/////////////////////////////////////////
/* The create function for the KD-Tree */
/////////////////////////////////////////
template<long D>
KDTree<D> *KDTreeCreate(const std::vector<IsotheticBBox<D> >& bbox);

/////////////////////////////////////////
/* The median function for the KD-Tree */
/////////////////////////////////////////
template<long D>
void MedianSelect(long , long , double *, long *, int);

template<long D>
void LocatePoint(const Point<D>& qp, 
        const KDTree<D>* kdtree,
        std::vector<long>& pfound);

template<long D>
void Intersect(const IsotheticBBox<D>& box, 
        const KDTree<D>* kdtree,
        std::vector<long>& pfound);


/****************************************************************************/
/* File           :MedianSelect.c                                           */
/*                                                                          */
/* Creation Date  :April 20, 1996                                           */
/*                                                                          */
/* Purpose        :MedianSelect takes an long `k', a double array 'arr',    */
/*                 and an long array `prm' of length `n'                    */
/*                 and reorders `prm' so that `arr[prm[k]]' is the k'th     */
/*                 largest element in the set `{arr[prm[i]], i=1,..,n}',    */
/*                 and we have in addition the partitioning properties:     */
/*                 i<k ==> arr[prm[i]] <= arr[prm[k]] and                   */
/*                 i>k ==> arr[prm[i]] >= arr[prm[k]].                      */
/*                                                                          */
/*                 Note that setting `k=(1+n)/2', this subroutine           */ 
/*                 computes the MEDIAN VALUE of `n' values in `arr',        */
/*                 and it does this in Order(N) time. This is thus a more   */
/*                 efficient way of computing the median than a heap sort   */
/*                 which is Order(NlogN).                                   */
/*                                                                          */
/****************************************************************************/

template<long D> void MedianSelect (long k, 
                                    long n, 
                                    std::vector<Point<D> >& arr, 
                                    long *prm, 
                                    long icut)
{
    long i, j, r, l, mid, ia;

    l = 0;
    r = n-1;

    while ( r-l > 1 ) {

        mid = (l+r)/2;
        SWAP(prm[mid], prm[l+1], long);

        if (arr[prm[l+1]][icut] > arr[prm[r]][icut]) { 
            SWAP(prm[l+1], prm[r], long);
        }

        if (arr[prm[l]][icut] > arr[prm[r]][icut]) {
            SWAP(prm[l], prm[r], long);
        }

        if (arr[prm[l+1]][icut] > arr[prm[l]][icut]) {
            SWAP(prm[l+1], prm[l], long);
        }

        i = l + 1;
        j = r;
        ia = prm[l];

        i++; 
        while ( arr[prm[i]][icut] < arr[ia][icut] ) i++;

        j--;
        while ( arr[prm[j]][icut] > arr[ia][icut] ) j--;

        while ( j >= i ) {
            SWAP(prm[i], prm[j], long);

            i++;
            while ( arr[prm[i]][icut] < arr[ia][icut] ) i++;

            j--;
            while ( arr[prm[j]][icut] > arr[ia][icut] ) j--;
        }

        prm[l] = prm[j];
        prm[j] = ia;

        if (j >= k) r = j - 1;
        if (j <= k) l = i;
    }

    if ( (r-l == 1) && (arr[prm[r]][icut] < arr[prm[l]][icut]) ) {
        SWAP(prm[l], prm[r], long);
    }

}




/****************************************************************************/
/* File           :kdtree.c                                                 */
/*                                                                          */
/* Creation Date  :April 20, 1996                                           */
/*                                                                          */
/*    MODIFIED: 15 March, 2001                                              */
/*                                                                          */
/* Purpose        :KDTREE takes the set of Safety Boxes and                 */
/*                 produces a k-D tree that is stored in the array LINKP.   */
/*                 Leaf nodes in LINKP each coincide with exactly one       */
/*                 Safety Box.  For each node in the k-D tree,              */
/*                 there is a corresponding Safety Box which is just        */
/*                 big enough to contain all the Safety Boxes ``under''     */
/*                 the node.                                                */
/*                                                                          */
/****************************************************************************/

template <long D>
KDTree<D> *KDTreeCreate(const std::vector<IsotheticBBox<D> >& sboxp)
{
    KDTree<D> *kdtree;
    std::vector<Point<D> > bbc;
    long *ipoly = NULL;
    long i;
    long node, nextp, itop;
    long imn, imx, icut, imd;
    long icrstack[100], imin[100], imax[100]; 
    long ict[100];
    Vector<D> dim;


    if( sboxp.empty()) std::abort(); 
    else {

        kdtree = new KDTree<D>;

        /* Obtain allocations for arrays */
        kdtree->num_entities = sboxp.size();
        ipoly = new long[sboxp.size()];
        kdtree->sbox = new IsotheticBBox<D>[2*sboxp.size()];
        kdtree->linkp = new long[2*sboxp.size()];
        bbc.resize(sboxp.size());

        for (i = 0; i < sboxp.size(); i++) {

            /* Compute the centers of the bounding boxes */
            bbc[i] = sboxp[i].center();

            /* Compute the root node of the k-D tree */
            kdtree->sbox[0].add(sboxp[i]);
        }

        /* If there is only one safety box, the root node is a leaf.
           (Our convention is to set the link corresponding to a leaf
           equal to the negative of the unique item contained in
           that leaf.)  If the root is a leaf, our work is done. */

        if ( sboxp.size() == 1 ) kdtree->linkp[0] = 0;
        else 
        {
            /* dimx, dimy, dimz are equal to the x, y, and z dimensions
               of the bounding box corresponding to the current node
               under consideration. */

            dim = kdtree->sbox[0].getMax() - kdtree->sbox[0].getMin();


            /* ipoly will contain a permutation of the integers
               {0,...,sboxp.size()-1}. This permutation will be altered as we
               create our balanced binary tree. nextp contains the
               address of the next available node to be filled. */

            nextp = 1;
            for (i = 0; i < sboxp.size(); i++) ipoly[i] = i;

            /* Use a stack to create the tree. Put the root node
               ``0'' on the top of the stack (icrstack). The array
               subset of ipoly (i.e., subset of boxes associated
               with this node) is recorded using the arrays
               imin and imax. */

            itop = 0;
            icrstack[itop] = 0;
            imin[itop] = 0;
            imax[itop] = sboxp.size()-1;

            /* Finally, we record in ict the appropriate ``cutting
               direction'' for bisecting the set of safety boxes
               corresponding to this node. This direction is either
               the x, y, or z directions, depending on which dimension
               of the bounding box is largest. */

            MaxComponent(dim,ict[itop]);

            /* Pop nodes off stack, create children nodes and put them
               on stack. Continue until k-D tree has been created. */

            while ( itop >= 0) { 

                /* Pop top node off stack. */

                node = icrstack[itop];

                /* Make this node point to next available node location (nextp).
                   This link represents the location of the FIRST CHILD
                   of the node. The adjacent location (nextp+1)
                   is implicitly taken to be the location of the SECOND
                   child of the node. */

                kdtree->linkp[node] = nextp;
                imn = imin[itop];
                imx = imax[itop];
                icut = ict[itop];
                itop--;

                /* Partition safety box subset associated with this node.
                   Using the appropriate cutting direction, use SELECT to
                   reorder (ipoly) so that the safety box with median bounding
                   box center coordinate is ipoly(imd), while the
                   boxes {ipoly[i], i<imd} have SMALLER (or equal)
                   bounding box coordinates, and the boxes with
                   {ipoly[i], i>imd} have GREATER (or equal) bounding box
                   coordinates. */

                imd = (imn+imx)/2;

                MedianSelect(imd-imn+1,imx-imn+1,bbc,&(ipoly[imn]),icut);

                /* If the first child's subset of safety boxes is a singleton,
                   the child is a leaf. Set the child's link to point to the
                   negative of the box number. Set the child's bounding
                   box to be equal to the safety box box. */

                if (imn == imd) {
                    kdtree->linkp[nextp]= -ipoly[imn]; 
                    kdtree->sbox[nextp] = sboxp[ipoly[imn]];
                    nextp = nextp+1;
                }
                else {

                    /* In this case, the subset of safety boxes corres to the
                       first child is more than one, and the child is
                       not a leaf. Compute the bounding box of this child to
                       be the smallest box containing all the associated safety'
                       boxes. */

                    kdtree->sbox[nextp] = sboxp[ipoly[imn]];

                    for (i=imn+1; i<=imd; i++) 
                        kdtree->sbox[nextp].add(sboxp[ipoly[i]]);

                    /* Put the first child onto the stack, noting the
                       associated triangle subset in imin and imax, and
                       putting the appropriate cutting direction in ict. */

                    dim = kdtree->sbox[nextp].getMax() 
                        - kdtree->sbox[nextp].getMin();

                    itop++;
                    icrstack[itop] = nextp;
                    imin[itop] = imn;
                    imax[itop] = imd;

                    MaxComponent(dim,ict[itop]);
                    nextp++;
                }

                /* If the second child's subset of safety boxes is a singleton,
                   the child is a leaf. Set the child's link to point to the
                   negative of the sbox number. Set the child's bounding
                   box to be equal to that of the safety box. */

                if ((imd+1) == imx) {
                    kdtree->linkp[nextp] = -ipoly[imx];
                    kdtree->sbox[nextp] = sboxp[ipoly[imx]];
                    nextp = nextp+1;
                }
                else {

                    /* In this case, the subset of boxes corresponding to the
                       second child is more than one safety box, and the child is
                       not a leaf. Compute the bounding box of this child to
                       be the smallest box containing all the associated safety
                       boxes. */

                    kdtree->sbox[nextp] = sboxp[ipoly[imd+1]];

                    for (i=imd+2; i<=imx; i++) 
                        kdtree->sbox[nextp].add(sboxp[ipoly[i]]);

                    /* Put the second child onto the stack, noting the
                       associated triangle subset in imin and imax, and
                       putting the appropriate cutting direction in ict. */

                    dim = kdtree->sbox[nextp].getMax() 
                        - kdtree->sbox[nextp].getMin();

                    itop++;
                    icrstack[itop] = nextp;
                    imin[itop] = imd+1;
                    imax[itop] = imx;

                    MaxComponent(dim,ict[itop]);
                    nextp++;
                }

            } /* End of the while loop */

            delete [] ipoly;
        }
    } 

    return kdtree;
}



// Return a list of BBox id's containing the query point
template<long D>
void LocatePoint(const Point<D>& qp, 
                 const KDTree<D>* kdtree, 
                 std::vector<long>& pfound)    
{
    pfound.clear();
    long itop, node, ind, j;
    long istack[100];

    /* If root node is a leaf, return leaf. */

    if (kdtree->linkp[0] <= 0) {
        if (kdtree->sbox[0].intersect(qp))
            pfound.push_back(-kdtree->linkp[0]);
    }
    else {        
        itop = 0;
        istack[itop] = 0;

        /* Traverse (relevant part of) k-D tree using stack. */

        while (itop >= 0) {

            /* pop node off of stack. */
            node = istack[itop]; itop--;

            ind = kdtree->linkp[node];

            /* check if either child of NODE is a leaf or should be
               put on stack. */
            for (j=0; j<=1; j++) { 
                if (kdtree->sbox[ind+j].intersect(qp)) {
                    /* If child is a leaf, add to the list. */
                    if (kdtree->linkp[ind+j] <= 0) { 	       
                        pfound.push_back(-kdtree->linkp[ind+j]);
                    }
                    else {
                        itop++;
                        istack[itop] = ind+j;
                    }
                }
            }
        }
    }
}

// Return a list (pfound) of BBox ids in the tree (kdtree) that overlap the 
// given BBox (box)
template<long D>
void Intersect(const IsotheticBBox<D>& box, 
               const KDTree<D>* kdtree, 
               std::vector<long>& pfound)    
{
    pfound.clear();
    long itop, node, ind, j;
    long istack[100];

    /* If root node is a leaf, return leaf. */
    if (kdtree->linkp[0] <= 0) {
        if (kdtree->sbox[0].intersect(box))
            pfound.push_back(-kdtree->linkp[0]);
    }
    else {        
        itop = 0;
        istack[itop] = 0;

        /* Traverse (relevant part of) k-D tree using stack. */
        while (itop >= 0) {

            /* pop node off of stack. */
            node = istack[itop]; itop--;

            ind = kdtree->linkp[node];

            /* check if either child of NODE is a leaf or should be
               put on stack. */
            for (j=0; j<=1; j++) { 
                if (kdtree->sbox[ind+j].intersect(box)) {

                    /* If child is a leaf, add to the list. */
                    if (kdtree->linkp[ind+j] <= 0) { 	       
                        pfound.push_back(-kdtree->linkp[ind+j]);
                    }
                    else {
                        itop++;
                        istack[itop] = ind+j;
                    }
                }
            }
        }
    }
}

#undef SWAP
}

#endif /* _KDTREE_H_ */

