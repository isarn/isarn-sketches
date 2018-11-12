/*
Copyright 2016-2018 Erik Erlandson

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

import org.scalatest._

import org.isarnproject.scalatest.matchers.seq._

class TDigestTest extends FlatSpec with Matchers {
  import org.apache.commons.math3.distribution.RealDistribution
  import org.apache.commons.math3.distribution.IntegerDistribution

  val seed = 235711L
  scala.util.Random.setSeed(seed)

  val ss = 100000
  val delta = 50.0 / 1000

  val maxD = 0.05
  val maxDI = 0.1

  def testTDvsDist(td: TDigest, dist: RealDistribution, stdv: Double): Boolean = {
    val xmin = td.clusters.keyMin.get
    val xmax = td.clusters.keyMax.get
    val step = (xmax - xmin) / 1000
    val d = (xmin to xmax by step).iterator
      .map(x => math.abs(td.cdf(x) - dist.cumulativeProbability(x))).max

    val dInv = (0.01 to 0.99 by 0.01).iterator
      .map(x => math.abs(td.cdfInverse(x) - dist.inverseCumulativeProbability(x))).max / stdv

    val pass = d <= maxD && dInv <= maxDI
    if (!pass) Console.err.println(s"testTDvsDist failure: d= $d  dInv= $dInv")
    pass
  }

  def testSamplingPDF(td: TDigest, dist: RealDistribution): Boolean = {
    val tdSamples = Array.fill(10000) { td.samplePDF }
    val distSamples = Array.fill(10000) { dist.sample }
    val kst = new org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest()
    val d = kst.kolmogorovSmirnovStatistic(tdSamples, distSamples)
    val pass = d <= maxD
    if (!pass) Console.err.println(s"testSamplingPDF failure: d= $d")
    pass
  }

  def testSamplingPMF(td: TDigest, dist: IntegerDistribution): Boolean = {
    td.nclusters should be <=(td.maxDiscrete)
    val tdSamples = Array.fill(10000) { td.samplePMF }
    val distSamples = Array.fill(10000) { dist.sample.toDouble }
    val kst = new org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest()
    val d = kst.kolmogorovSmirnovStatistic(tdSamples, distSamples)
    val pass = d <= maxD
    if (!pass) Console.err.println(s"testSamplingPDF failure: d= $d")
    pass
  }

  def testDistribution(dist: RealDistribution, stdv: Double): Boolean = {
    dist.reseedRandomGenerator(seed)

    val td = TDigest.sketch(Iterator.fill(ss) { dist.sample }, delta = delta)

    testTDvsDist(td, dist, stdv) && testSamplingPDF(td, dist)
  }

  def testMonotoneCDF(dist: RealDistribution): Boolean = {
    dist.reseedRandomGenerator(seed)
    val td = TDigest.sketch(Iterator.fill(ss) { dist.sample }, delta = delta)
    val (xmin, xmax) = (td.clusters.keyMin.get, td.clusters.keyMax.get)
    val step = (xmax - xmin) / 100000
    val t = (xmin to xmax by step).iterator.map(x => td.cdf(x)).sliding(2).map(w => w(1) - w(0)).min
    val pass = t >= 0.0
    if (!pass) Console.err.println(s"testMonotoneCDF failure: t= $t")
    pass
  }

  def testMonotoneCDFI(dist: RealDistribution): Boolean = {
    dist.reseedRandomGenerator(seed)
    val td = TDigest.sketch(Iterator.fill(ss) { dist.sample }, delta = delta)
    val (xmin, xmax) = (0.0, 1.0)
    val step = (xmax - xmin) / 100000
    val t = (xmin to xmax by step).iterator.map(q => td.cdfInverse(q)).sliding(2).map(w => w(1) - w(0)).min
    val pass = t >= 0.0
    if (!pass) Console.err.println(s"testMonotoneCDFI failure: t= $t")
    pass
  }

  def testMonotone(dist: RealDistribution): Boolean = {
    testMonotoneCDF(dist) && testMonotoneCDFI(dist)
  }

  it should "sketch a uniform distribution" in {
    import org.apache.commons.math3.distribution.UniformRealDistribution
    val dist = new UniformRealDistribution()
    testDistribution(dist, math.sqrt(dist.getNumericalVariance())) should be (true)
  }

  it should "sketch a normal distribution" in {
    import org.apache.commons.math3.distribution.NormalDistribution
    val dist = new NormalDistribution()
    testDistribution(dist, math.sqrt(dist.getNumericalVariance())) should be (true)
  }

  it should "sketch an exponential distribution" in {
    import org.apache.commons.math3.distribution.ExponentialDistribution
    val dist = new ExponentialDistribution(1.0)
    testDistribution(dist, math.sqrt(dist.getNumericalVariance())) should be (true)
  }

  it should "aggregate with another t-digest using ++" in {
    import org.apache.commons.math3.distribution.NormalDistribution
    val dist = new NormalDistribution()
    dist.reseedRandomGenerator(seed)

    val td1 = TDigest.sketch(Iterator.fill(ss) { dist.sample }, delta = delta)
    val td2 = TDigest.sketch(Iterator.fill(ss) { dist.sample }, delta = delta)

    testTDvsDist(td1 ++ td2, dist, math.sqrt(dist.getNumericalVariance())) should be (true)
  }

  it should "respect monotonic cdf and inverse" in {
    import org.apache.commons.math3.distribution.ExponentialDistribution
    import org.apache.commons.math3.distribution.NormalDistribution
    import org.apache.commons.math3.distribution.UniformRealDistribution

    testMonotone(new UniformRealDistribution()) should be (true)
    testMonotone(new ExponentialDistribution(1.0)) should be (true)
    testMonotone(new NormalDistribution(0.0, 0.1)) should be (true)
  }

  it should "respect maxDiscrete parameter" in {
    import org.apache.commons.math3.distribution.GeometricDistribution
    val gd = new GeometricDistribution(0.33)
    val data = gd.sample(1000000)
    val dataUniq = data.distinct.sorted
    val kt = dataUniq.map(_.toDouble).toSet
    val td = TDigest.sketch(data, maxDiscrete = 50)
    val clust = td.clusters
    clust.keys.toSet should be (kt)
    val D = clust.keys.map { x => td.cdfDiscrete(x) }
      .zip(dataUniq.map { k => gd.cumulativeProbability(k) })
      .map { case (p1, p2) => math.abs(p1 - p2) }
      .max
    (D <= 0.01) should be (true)
    testSamplingPMF(td, gd) should be (true)
  }

  it should "respect maxDiscrete parameter over ++" in {
    import org.apache.commons.math3.distribution.GeometricDistribution
    val gd = new GeometricDistribution(0.33)
    val tdvec = Vector.fill(10) { TDigest.sketch(gd.sample(100000), maxDiscrete = 50) }
    val td = tdvec.reduce(_ ++ _)
    val clust = td.clusters
    clust.keys.map(_.toInt).map(_.toDouble) should beEqSeq(clust.keys)
    val D = clust.keys.map { x => td.cdfDiscrete(x) }
      .zip(clust.keys.map(_.toInt).map { k => gd.cumulativeProbability(k) })
      .map { case (p1, p2) => math.abs(p1 - p2) }
      .max
    (D <= 0.01) should be (true)
    testSamplingPMF(td, gd) should be (true)
  }

  it should "serialize and deserialize" in {
    import org.apache.commons.math3.distribution.NormalDistribution

    import org.isarnproject.scalatest.serde.roundTripSerDe

    val dist = new NormalDistribution()
    dist.reseedRandomGenerator(seed)

    val tdo = TDigest.sketch(Iterator.fill(ss) { dist.sample }, delta = delta)

    val tdi = roundTripSerDe(tdo)

    (tdi == tdo) should be (true)

    testTDvsDist(tdi, dist, math.sqrt(dist.getNumericalVariance())) should be (true)
  }
}
