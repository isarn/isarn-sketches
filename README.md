# isarn-sketches
Sketching data structures

### How to use in your project

#### sbt
``` scala
resolvers += "isarn project" at "https://dl.bintray.com/isarn/maven/"

libraryDependencies += "org.isarnproject" %% "isarn-sketches" % "0.0.1"
```

#### maven
``` xml
<dependency> 
  <groupId>org.isarnproject</groupId>
  <artifactId>isarn-sketches_2.10</artifactId> 
  <version>0.0.1</version> 
  <type>pom</type> 
</dependency>
```

### t-digest
``` scala
scala> import org.isarnproject.sketches.TDigest
import org.isarnproject.sketches.TDigest

scala> val data = Vector.fill(10000) { scala.util.Random.nextGaussian() }
data: scala.collection.immutable.Vector[Double] = Vector(1.6046163970051968, 0.44151418924289004, ...

scala> val sketch = TDigest.sketch(data)
sketch: org.isarnproject.sketches.TDigest = TDigest(0.5,70,TDigestMap(-3.6035923746624587 -> (1.0, 1.0), ...

scala> sketch.cdf(0)
res0: Double = 0.4984362744530557

scala> sketch.cdfInverse(0.5)
res1: Double = 0.0038481195948969205
```
