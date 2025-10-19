/*
 * =============================================================================
 *
 *       Filename:  rbtree.h
 *
 *    Description:  rbtree(Red-Black tree) implementation adapted from linux
 *                  kernel thus can be used in userspace c program.
 *
 *        Created:  09/02/2012 11:36:11 PM
 *
 *         Author:  Fu Haiping (forhappy), haipingf@gmail.com
 *        Company:  ICT ( Institute Of Computing Technology, CAS )
 *
 * =============================================================================
 */

/*
  Red Black Trees
  (C) 1999  Andrea Arcangeli <andrea@suse.de>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  linux/include/linux/rbtree.h

  To use rbtrees you'll have to implement your own insert and search cores.
  This will avoid us to use callbacks and to drop drammatically performances.
  I know it's not the cleaner way,  but in C (not in C++) to get
  performances and genericity...

  Some example of insert and search follows here. The search is a plain
  normal search over an ordered tree. The insert instead must be implemented
  in two steps: First, the code must insert the element in order as a red leaf
  in the tree, and then the support library function rb_insert_color() must
  be called. Such function will do the not trivial work to rebalance the
  rbtree, if necessary.
*/

#ifndef _LINUX_RBTREE_H
#define _LINUX_RBTREE_H

#if defined(container_of)
#undef container_of
#define container_of(ptr, type, member) ((type *)(((char *)(ptr)) - ((char *)(&((type *)0)->member))))
#else
#define container_of(ptr, type, member) ((type *)(((char *)(ptr)) - ((char *)(&((type *)0)->member))))
#endif

#undef NULL
#if defined(__cplusplus)
#define NULL 0
#else
#define NULL ((void *)0)
#endif

struct rb_node
{
    unsigned long rb_parent_color;
#define RB_RED   0
#define RB_BLACK 1
    struct rb_node *rb_right;
    struct rb_node *rb_left;
} __attribute__((aligned(sizeof(long))));
/* The alignment might seem pointless, but allegedly CRIS needs it */

struct rb_root
{
    struct rb_node *rb_node;
};

typedef struct rb_root *p_rb_root;
typedef void (*rb_augment_f)(struct rb_node *node, void *data);

#endif /* _LINUX_RBTREE_H */
