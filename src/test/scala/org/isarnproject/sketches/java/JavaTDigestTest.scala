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

package org.isarnproject.sketches.java

import org.scalatest._

import org.isarnproject.scalatest.matchers.seq._

class JavaTDigestTest extends FlatSpec with Matchers {
  import org.apache.commons.math3.distribution.RealDistribution
  import org.apache.commons.math3.distribution.IntegerDistribution

  val seed = 235711L
  scala.util.Random.setSeed(seed)

  val ss = 100000
  val delta = 50.0 / 1000

  val maxD = 0.05
  val maxDI = 0.1

  def testTDvsDist(td: TDigest, dist: RealDistribution, stdv: Double): Boolean = {
    val xmin = td.cent(0)
    val xmax = td.cent(td.nclusters - 1)
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

    val td = TDigest.sketch(Array.fill(ss) { dist.sample }, delta)

    testTDvsDist(td, dist, stdv) && testSamplingPDF(td, dist)
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

  it should "aggregate with another t-digest using merge method" in {
    import org.apache.commons.math3.distribution.NormalDistribution
    val dist = new NormalDistribution()
    dist.reseedRandomGenerator(seed)

    val td1 = TDigest.sketch(Array.fill(ss) { dist.sample }, delta)
    val td2 = TDigest.sketch(Array.fill(ss) { dist.sample }, delta)

    testTDvsDist(TDigest.merge(td1, td2), dist, math.sqrt(dist.getNumericalVariance())) should be (true)
  }
}
