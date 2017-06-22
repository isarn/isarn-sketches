/*
Copyright 2016-2017 Erik Erlandson

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

package org.isarnproject.sketches.tdmap

import math.Numeric

import scala.collection.SortedMap

import org.isarnproject.algebraAPI.{ MonoidAPI => Monoid, AggregatorAPI => Aggregator }

import org.isarnproject.collections.mixmaps.increment._
import org.isarnproject.collections.mixmaps.prefixsum._
import org.isarnproject.collections.mixmaps.nearest._

object tree {
  import org.isarnproject.collections.mixmaps.redblack.tree._
  import org.isarnproject.collections.mixmaps.ordered._
  import org.isarnproject.collections.mixmaps.ordered.tree.DataMap
  import org.isarnproject.collections.mixmaps.increment.tree._
  import org.isarnproject.collections.mixmaps.prefixsum.tree._
  import org.isarnproject.collections.mixmaps.nearest.tree._

  /** Base class of tree node for a TDigestMap */
  trait NodeTD extends NodePS[Double, Double, Double]
    with NodeInc[Double, Double] with NodeNearMap[Double, Double] {

    /**
     * Obtain a "mass cover": two adjacent nodes in the tree such that the cumulative mass
     * of the left node is <= (m) and the cumulative mass of the right node is > (m)
     * @param m The target mass to cover between two adjacent nodes
     * @return a Cover instance with the left and right covering tree nodes.  If the (m)
     * was < the mass of the left-most tree node, the left cover value will be None.  Similarly
     * if the mass was >= the cumulative mass of the right-most node (equivalent to sum of all
     * node masses in the tree), then the right cover value will be None.
     */
    final def mCover(m: Double) = mcov(m, 0.0, Cover[INodeTD](None, None))

    private[tree] def mcov(m: Double, psum: Double, cov: Cover[INodeTD]): Cover[INodeTD]

    // Find the cluster whose prefix sum is the least upper bound of mass 'm'
    // domain specific to t-digest algorithms
    private[sketches] final def keyPFSLUB(m: Double) = this match {
      case _: LNodeTD => Double.NaN
      case _ if (m < 0.0 || m > this.pfs) => Double.NaN
      case _ if (m == 0.0) => this.nodeMin.get.asInstanceOf[INodeTD].data.key
      case _ => this.kpl(m, 0.0)
    }

    // recursive implementation of keyPFSLUB
    private[tree] def kpl(m: Double, psum: Double): Double

    // obtains the nearest cluster to 'x'.  Returns the cluster (location, mass, prefix-sum)
    private[sketches] final def nearTD(x: Double): (Double, Double, Double) = ntd(x, 0.0)

    // recursive implementation for nearTD
    private[tree] def ntd(x: Double, psum: Double): (Double, Double, Double)

    // recursive implementation of 'update' method
    private[tdmap] def upd(x0: Double, x: Double, m: Double): Node[Double]
  }

  trait LNodeTD extends NodeTD
    with LNodePS[Double, Double, Double] with LNodeInc[Double, Double]
    with LNodeNearMap[Double, Double] {
    final def mcov(m: Double, psum: Double, cov: Cover[INodeTD]) = cov
    final def kpl(m: Double, psum: Double) = Double.NaN
    final def ntd(x: Double, psum: Double) = (Double.NaN, Double.NaN, Double.NaN)
    final def upd(x0: Double, x: Double, m: Double) =
      throw new Exception("If this exception threw, there is a bug in this code")
  }

  trait INodeTD extends NodeTD
    with INodePS[Double, Double, Double] with INodeInc[Double, Double]
    with INodeNearMap[Double, Double] {
    val lsub: NodeTD
    val rsub: NodeTD

    final def mcov(m: Double, psum: Double, cov: Cover[INodeTD]) = {
      if (m < psum + lsub.pfs) {
        lsub match {
          case n: INodeTD =>
            lsub.mcov(m, psum, cov.copy(r = Some(n.nodeMax.get.asInstanceOf[INodeTD])))
          case _ => cov.copy(r = Some(this))
        }
      } else {
        val t = psum + lsub.pfs + data.value
        if (m >= t) rsub.mcov(m, t, cov.copy(l = Some(this)))
        else {
          lsub match {
            case n: INodeTD => Cover(Some(lsub.nodeMax.get.asInstanceOf[INodeTD]), Some(this))
            case _ => cov.copy(r = Some(this))
          }
        }
      }
    }

    final def kpl(m: Double, psum: Double) = {
      val lb = psum + lsub.pfs
      val ub = lb + data.value
      if (m > ub) {
        rsub.kpl(m, ub)
      } else if (m > lb) {
        data.key
      } else {
        lsub.kpl(m, psum)
      }
    }

    final def ntd(x: Double, psum: Double) = {
      if (x < data.key) {
        lsub match {
          case ls: INodeTD => {
            if (x <= ls.kmax) ls.ntd(x, psum)
            else {
              val (dk, ldk) = (math.abs(x - data.key), math.abs(x - ls.kmax))
              if (dk <= ldk) (data.key, data.value, psum + lsub.pfs)
              else {
                val n = ls.node(ls.kmax).get.asInstanceOf[INodeTD]
                (n.data.key, n.data.value, psum + lsub.pfs - n.data.value)
              }
            }
          }
          case _ => (data.key, data.value, psum + lsub.pfs)
        }
      } else if (x > data.key) {
        rsub match {
          case rs: INodeTD => {
            if (x >= rs.kmin) rs.ntd(x, psum + lsub.pfs + data.value)
            else {
              val (dk, rdk) = (math.abs(x - data.key), math.abs(x - rs.kmin))
              if (dk <= rdk) (data.key, data.value, psum + lsub.pfs)
              else {
                val n = rs.node(rs.kmin).get.asInstanceOf[INodeTD]
                (n.data.key, n.data.value, psum + lsub.pfs + data.value)
              }
            }
          }
          case _ => (data.key, data.value, psum + lsub.pfs)
        }
      } else (data.key, data.value, psum + lsub.pfs)
    }

    final def upd(x0: Double, x: Double, m: Double) =
      if (color == R) {
        if (x0 < data.key) rNode(data, lsub.upd(x0, x, m), rsub)
        else if (x0 > data.key) rNode(data, lsub, rsub.upd(x0, x, m))
        else {
          val d = new DataMap[Double, Double] {
            val key = x
            val value = m
          }
          rNode(d, lsub, rsub)
        }
      } else {
        // We know we are directly replacing a node, so no need to call balance()
        // in the case of black nodes. This is quite a bit faster. \o/
        if (x0 < data.key) bNode(data, lsub.upd(x0, x, m), rsub)
        else if (x0 > data.key) bNode(data, lsub, rsub.upd(x0, x, m))
        else {
          val d = new DataMap[Double, Double] {
            val key = x
            val value = m
          }
          bNode(d, lsub, rsub)
        }
      }
  }
}

import tree._

object infra {
  import org.isarnproject.collections.mixmaps.redblack.tree._
  import org.isarnproject.collections.mixmaps.ordered.tree.DataMap

  object tdmapMonoid extends Monoid[Double] {
    def empty = 0.0
    def combine(x: Double, y: Double) = x + y
    def combineAll(as: TraversableOnce[Double]) = as.fold(0.0)(_ + _)
    def combineAllOption(as: TraversableOnce[Double]) =
      if (as.isEmpty) None else Some(combineAll(as))
  }

  object tdmapAggregator extends Aggregator[Double, Double] {
    def monoid = tdmapMonoid
    def lff = (m: Double, d: Double) => m + d
    def mf = (d: Double) => d
    def aggregate(as: TraversableOnce[Double]) = as.fold(0.0)(_ + _)
  }

  /** Dependency injection class for TDigestMap */
  class Inject extends Serializable {
    // Typeclasses corresponding to "regular real numbers":
    val keyOrdering = implicitly[Numeric[Double]]

    val valueMonoid = tdmapMonoid
    
    val prefixAggregator = tdmapAggregator

    def iNode(clr: Color, dat: Data[Double], ls: Node[Double], rs: Node[Double]) =
      new Inject with INodeTD with TDigestMap {
        // INode
        val color = clr
        val lsub = ls.asInstanceOf[NodeTD]
        val rsub = rs.asInstanceOf[NodeTD]
        val data = dat.asInstanceOf[DataMap[Double, Double]]
        // INodePS
        val prefix = prefixAggregator.lff(
          prefixAggregator.monoid.combine(lsub.pfs, rsub.pfs), data.value)
        // INodeNear
        val kmin = lsub match {
          case n: INodeTD => n.kmin
          case _ => data.key
        }
        val kmax = rsub match {
          case n: INodeTD => n.kmax
          case _ => data.key
        }
      }
  }

}

import infra._

/**
 * The tree-backed map object a TDigest uses to store and update its clusters.  TDigestMap
 * inherits functionality for value increment, prefix-sum and nearest-neighbor queries.
 */
sealed trait TDigestMap extends SortedMap[Double, Double] with NodeTD
  with IncrementMapLike[Double, Double, INodeTD, TDigestMap]
  with PrefixSumMapLike[Double, Double, Double, INodeTD, TDigestMap]
  with NearestMapLike[Double, Double, INodeTD, TDigestMap] {

  override def empty = TDigestMap.empty

  private def m1m2(c1: Double, tm1: Double, c2: Double, tm2: Double) = {
    val s = this.prefixSum(c1, open = true)
    val d1 = if (c1 == this.keyMin.get) 0.0 else tm1 / 2.0
    val m1 = s + d1
    val m2 = m1 + (tm1 - d1) + (if (c2 == this.keyMax.get) tm2 else tm2 / 2.0)
    (m1, m2)
  }

  // This updates an existing cluster with a new location and mass.  It does this
  // efficiently by taking advantage of the knowledge that (a) this kind of update
  // never changes the key ordering, and therefore that (b) this operation can
  // always directly replace an existing node, without otherwise changing the topology
  // of the tree.  Clearly, this is a domain-dependent method, and not exposed to the
  // public API
  private[sketches] def update(x0: Double, x: Double, m: Double): TDigestMap =
    this.upd(x0, x, m).asInstanceOf[TDigestMap]

  /** Compute the CDF for a value, using piece-wise linear between clusters */
  def cdf[N](xx: N)(implicit num: Numeric[N]) = {
    val x = num.toDouble(xx)
    this.coverR(x) match {
      case Cover(Some((c1, tm1)), Some((c2, tm2))) => {
        val (m1, m2) = m1m2(c1, tm1, c2, tm2)
        (m1 + (x - c1) * (m2 - m1) / (c2 - c1)) / this.sum
      }
      case Cover(Some(_), None) => 1.0
      case _ => 0.0
    }
  }

  def cdfDiscrete[N](xx: N)(implicit num: Numeric[N]) = {
    if (this.isEmpty) 0.0 else {
      val x = num.toDouble(xx)
      this.prefixSum(x) / this.sum
    }
  }

  def cdfDiscreteInverse[N](qq: N)(implicit num: Numeric[N]) = {
    val q = num.toDouble(qq)
    keyPFSLUB(q * this.sum)
  }

  /**
   * Compute the inverse-CDF from a given quantile on interval [0, 1], using piecewise linear
   * interpolation between clusters
   */
  def cdfInverse[N](qq: N)(implicit num: Numeric[N]) = {
    def cdfI(m: Double, c1: Double, tm1: Double, c2: Double, tm2: Double) = {
      val (m1, m2) = m1m2(c1, tm1, c2, tm2)
      c1 + (m - m1) * (c2 - c1) / (m2 - m1)
    }

    val q = num.toDouble(qq)
    if (q < 0.0 || q > 1.0) Double.NaN
    else {
      val m = q * this.sum
      this.mCover(m).map(n => (n.data.key, n.data.value)) match {
        case Cover(Some((c1, tm1)), Some((c2, tm2))) => cdfI(m, c1, tm1, c2, tm2)
        case Cover(None, Some((c, _))) => this.coverR(c) match {
          case Cover(Some((c1, tm1)), Some((c2, tm2))) => cdfI(m, c1, tm1, c2, tm2)
          case _ => Double.NaN
        }
        case Cover(Some((c, _)), None) => c
        case _ => Double.NaN
      }
    }
  }

  override def toString =
    "TDigestMap(" +
      iterator.zip(prefixSumsIterator())
      .map(x => s"${x._1._1} -> (${x._1._2}, ${x._2})").mkString(", ") +
      ")"
}

/** factory functions for TDigestMap */
object TDigestMap {
  /** Obtain an empty TDigestMap instance */
  def empty: TDigestMap = new Inject with LNodeTD with TDigestMap
}
