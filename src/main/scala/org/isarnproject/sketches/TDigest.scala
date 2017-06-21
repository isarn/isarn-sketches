/*
Copyright 2016 Erik Erlandson

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

package org.isarnproject.sketches

import scala.util.Random

import tdmap.TDigestMap

/**
 * A t-digest sketch of sampled numeric data, as described in:
 * Computing Extremely Accurate Quantiles Using t-Digests,
 * Ted Dunning and Otmar Ertl,
 * https://github.com/tdunning/t-digest/blob/master/docs/t-digest-paper/histo.pdf
 *
 * {{{
 * import org.isarnproject.sketches.TDigest
 * val data = Vector.fill(10000) { scala.util.Random.nextGaussian() }
 * // sketch of some Gaussian data
 * val sketch = TDigest.sketch(data)
 * // the cumulative distribution function of the sketch; cdf(x) at x = 0
 * val cdf = sketch.cdf(0.0)
 * // inverse of the CDF, evaluated at q = 0.5
 * val cdfi = sketch.cdfInverse(0.5)
 * }}}
 */
case class TDigest(
  delta: Double,
  maxDiscrete: Int,
  nclusters: Int,
  clusters: TDigestMap) extends Serializable {

  // re-cluster when number of clusters exceeds this threshold
  private val R = (TDigest.K / delta).toInt

  private case class Cluster(centroid: Double, mass: Double, massUB: Double)

  /**
   * Returns a new t-digest with value x included in its sketch; td + x is equivalent to
   * td + (x, 1).
   * @param x The numeric data value to include in the sketch
   * @return the updated sketch
   */
  def +[N](x: N)(implicit num: Numeric[N]): TDigest = this.plus(num.toDouble(x), 1.0)

  /**
   * Returns a new t-digest with new pair (x, w) included in its sketch.
   * @param xw A pair (x, w) where x is the numeric value and w is its weight
   * @return the updated sketch
   * @note This implements 'algorithm 1' from:
   * Computing Extremely Accurate Quantiles Using t-Digests,
   * Ted Dunning and Otmar Ertl,
   * https://github.com/tdunning/t-digest/blob/master/docs/t-digest-paper/histo.pdf
   */
  def +[N1, N2](xw: (N1, N2))(implicit num1: Numeric[N1], num2: Numeric[N2]): TDigest =
    this.plus(num1.toDouble(xw._1), num2.toDouble(xw._2))

  private def plus(x: Double, w: Double): TDigest = {
    if (nclusters <= maxDiscrete) {
      val ncNew = nclusters + (if (clusters.contains(x)) 0 else 1)
      TDigest(delta, maxDiscrete, ncNew, clusters.increment(x, w))
    } else {
      val s = this.update(x, w)
      if (s.nclusters <= R) s
      else {
        // too many clusters: attempt to compress it by re-clustering
        val ds = TDigest.shuffle(s.clusters.toVector)
        ds.foldLeft(TDigest.empty(delta, maxDiscrete)) { case (d, (x, w)) => d.update(x, w) }
      }
    }
  }

  /**
   * Add this digest to another
   * @param that The right-hand t-digest operand
   * @return the result of combining left and right digests
   */
  def ++(that: TDigest): TDigest = TDigest.combine(this, that, this.delta, this.maxDiscrete)

  // This is most of 'algorithm 1', except for re-clustering which is factored out to avoid
  // recursive calls during a reclustering phase
  private def update(x: Double, w: Double) = {
    require(w > 0.0, "data weight must be > 0")

    if (clusters.isEmpty) {
      // our map is empty, so insert this pair as the first cluster
      TDigest(delta, maxDiscrete, nclusters + 1, clusters + (x -> w))
    } else {
      // Get the current cluster nearest to incoming (x)
      val (c, m, psum) = clusters.nearTD(x)
      if (x == c) {
        // data landed on an existing cluster: increment that cluster's mass directly
        TDigest(delta, maxDiscrete, nclusters, clusters.increment(c, w))
      } else {
        val M = clusters.sum
        val q = (psum + m / 2.0) / M
        val ub = M * delta * q * (1.0 - q)

        val dm = math.min(w, math.max(0.0, ub - m))
        val rm = w - dm

        val tClust = if (dm > 0.0) {
          val nm = m + dm
          val dc = dm * (x - c) / nm
          (clusters - c) + ((c + dc) -> nm)
        } else clusters

        val uClust = if (rm > 0.0) tClust + (x -> rm) else tClust

        // return the updated t-digest
        TDigest(delta, maxDiscrete, nclusters + (if (rm > 0.0) 1 else 0), uClust)
      }
    }
  }

  /**
   * Compute a cumulative probability (CDF) for a numeric value, from the estimated probability
   * distribution represented by this t-digest sketch.
   * @param x a numeric value
   * @return the cumulative probability that a random sample from the distribution is <= x
   */
  def cdf[N](x: N)(implicit num: Numeric[N]): Double = clusters.cdf(x)

  /**
   * Compute the inverse cumulative probability (inverse-CDF) for a quantile value, from the
   * estimated probability distribution represented by this t-digest sketch.
   * @param q a quantile value.  The value of q is expected to be on interval [0, 1]
   * @return the value x such that cdf(x) = q
   */
  def cdfInverse[N](q: N)(implicit num: Numeric[N]): Double = clusters.cdfInverse(q)

  /**
   * Compute a cumulative probability (CDF) for a numeric value, from the estimated probability
   * distribution represented by this t-digest sketch, assuming sketch is "discrete"
   * (e.g. if number of clusters <= maxDiscrete setting)
   * @param x a numeric value
   * @return the cumulative probability that a random sample from the distribution is <= x
   */
  def cdfDiscrete[N](x: N)(implicit num: Numeric[N]): Double =
    clusters.cdfDiscrete(x)

  /**
   * Compute the inverse cumulative probability (inverse-CDF) for a quantile value, from the
   * estimated probability distribution represented by this t-digest sketch,
   * assuming the sketch is "discrete" (e.g. if number of clusters <= maxDiscrete setting)
   * @param q a quantile value.  The value of q is expected to be on interval [0, 1]
   * @return the smallest value x such that q <= cdf(x)
   */
  def cdfDiscreteInverse[N](q: N)(implicit num: Numeric[N]): Double =
    clusters.cdfDiscreteInverse(q)

  /**
   * Perform a random sampling from the distribution as sketched by this t-digest, in
   * "probability density" mode.
   * @return A random number sampled from the sketched distribution
   * @note uses the inverse transform sampling method
   */
  def samplePDF: Double = clusters.cdfInverse(Random.nextDouble)

  /**
   * Perform a random sampling from the distribution as sketched by this t-digest, in
   * "probability mass" (i.e. discrete) mode.
   * @return A random number sampled from the sketched distribution
   * @note uses the inverse transform sampling method
   */
  def samplePMF: Double = clusters.cdfDiscreteInverse(Random.nextDouble)

  /**
   * Perform a random sampling from the distribution as sketched by this t-digest,
   * using "discrete" (PMF) mode if the number of clusters <= maxDiscrete setting,
   * and "density" (PDF) mode otherwise.
   * @return A random number sampled from the sketched distribution
   * @note uses the inverse transform sampling method
   */
  def sample: Double = if (nclusters <= maxDiscrete) samplePMF else samplePDF
}

/** Factory functions for TDigest */
object TDigest {
  import scala.language.higherKinds
  import scala.collection.SeqLike
  import scala.collection.generic.CanBuildFrom

  /**
   * Default value for a t-digest delta parameter.  The number of clusters varies, roughly, as
   * about (50/delta), when data are presented in random order
   * (it may grow larger if data are not presented randomly).  The default corresponds to
   * an expected number of clusters of about 100.
   */
  val deltaDefault: Double = (50.0 / 100.0) // delta * E[clusters] ~ 50

  /**
   * The t-digest algorithm will re-cluster itself whenever its number of clusters exceeds
   * (K/delta).  This value is set such that the threshold is about 10x the heuristically
   * expected number of clusters for the user-specified delta value.  Generally the number of
   * clusters will only trigger the corresponding re-clustering threshold when data are being
   * presented in a non-random order.
   */
  val K: Double = 10.0 * 50.0

  /**
   * Obtain an empty t-digest
   * @param delta a sketch resolution parameter.
   * @note Smaller values of delta yield sketches with more clusters, and higher resolution
   * @note The expected number of clusters will vary (roughly) as (50/delta)
   */
  def empty(delta: Double = deltaDefault, maxDiscrete: Int = 0): TDigest = {
    require(delta > 0.0, s"delta was not > 0")
    require(maxDiscrete >= 0, s"maxDiscrete was not >= 0")
    TDigest(delta, maxDiscrete, 0, TDigestMap.empty)
  }

  /**
   * Sketch some data with a t-digest
   * @param data The data elements to sketch
   * @param delta The sketch resolution parameter.
   * @return A t-digest sketch of the input data
   * @note Smaller values of delta yield sketches with more clusters, and higher resolution
   * @note The expected number of clusters will vary (roughly) as (50/delta)
   */
  def sketch[N](
    data: TraversableOnce[N],
    delta: Double = deltaDefault,
    maxDiscrete: Int = 0)(implicit num: Numeric[N]): TDigest = {
    require(delta > 0.0, s"delta was not > 0")
    require(maxDiscrete >= 0, s"maxDiscrete was not >= 0")
    val td = data.foldLeft(empty(delta, maxDiscrete))((c, e) => c + ((e, 1)))
    TDigest.shuffle(td.clusters.toVector).foldLeft(empty(delta, maxDiscrete))((c, e) => c + e)
  }

  /**
   * Combine two t-digests to yield a new digest
   * @param ltd the left-hand t-digest operand
   * @param rtd the right hand t-digest
   * @return the sum of left and right digests, defined as their aggregation
   * @note This operation satisfies a Semigroup law, with the caveat
   * that it is only "statistically" associative: d1++(d2++d3) will be statistically
   * similar to (d1++d2)++d3, but rarely identical.
   */
  def combine(ltd: TDigest, rtd: TDigest,
      delta: Double = deltaDefault,
      maxDiscrete: Int = 0): TDigest = {
    if (ltd.nclusters <= 1 && rtd.nclusters > 1) combine(rtd, ltd, delta, maxDiscrete)
    else if (rtd.nclusters == 0) ltd
    else if (rtd.nclusters == 1) {
      // handle the singleton RHS case specially to prevent quadratic catastrophe when
      // it is being used in the Aggregator use case
      val d = rtd.clusters.asInstanceOf[tdmap.tree.INodeTD].data
      ltd + ((d.key, d.value))
    } else {
      // insert clusters from largest to smallest
      (ltd.clusters.toVector ++ rtd.clusters.toVector).sortWith((a, b) => a._2 > b._2)
        .foldLeft(empty(delta, maxDiscrete))((d, e) => d + e)
    }
  }

  // Shuffle a sequence in a referentially-transparent way: pseudo-randomly, but with a random
  // seed that is a function of the sequence argument.
  private[sketches] def shuffle[E, S[X] <: SeqLike[X, S[X]]](
    seq: S[E])(implicit cbf: CanBuildFrom[S[E], E, S[E]]) =
    if (seq.length <= 1) seq
    else {
      val seed = scala.util.hashing.MurmurHash3.productHash((seq(0), seq(1), seq.length))
      (new scala.util.Random(seed)).shuffle(seq)
    }
}
