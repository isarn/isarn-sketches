# isarn-sketches
Sketching data structures

### API documentation
https://isarn.github.io/isarn-sketches/scala/api/
https://isarn.github.io/isarn-sketches/java/api/

### Compatibility
isarn-sketches can operate with [Algebird](https://twitter.github.io/algebird/) via the
[isarn-sketches-algebird-api](https://github.com/isarn/isarn-sketches-algebird-api)

isarn-sketches can also operate with [Apache Spark](https://github.com/apache/spark) via the [isarn-sketches-spark](https://github.com/isarn/isarn-sketches-spark) library

### How to use in your project

``` scala
libraryDependencies += "org.isarnproject" %% "isarn-sketches" % "0.1.2"
```

### t-digest
``` scala
scala> import org.isarnproject.sketches.TDigest
import org.isarnproject.sketches.TDigest

scala> val data = Vector.fill(10000) { scala.util.Random.nextGaussian() }
data: scala.collection.immutable.Vector[Double] = Vector(1.6046163970051968, 0.44151418924289004, ...

scala> val sketch = TDigest.sketch(data)
sketch: org.isarnproject.sketches.TDigest = TDigest(0.5,0,74,TDigestMap(-3.819069044174932 -> (1.0, 1.0), ...

scala> sketch.cdf(0)
res0: Double = 0.4984362744530557

scala> sketch.cdfInverse(0.5)
res1: Double = 0.0038481195948969205
```

#### t-digest resources
* Original paper: [Computing Extremely Accurate Quantiles Using t-Digests](https://github.com/tdunning/t-digest/blob/master/docs/t-digest-paper/histo.pdf)
* Video Talk: [Sketching Data with T-Digest In Apache Spark](https://youtu.be/ETUYhEZRtWE)
