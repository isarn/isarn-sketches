resolvers ++= Seq(
  "jgit-repo".at("https://download.eclipse.org/jgit/maven"),
  //"sonatype-releases" at "https://oss.sonatype.org/content/repositories/releases/",
  //Resolver.url("bintray-sbt-plugin-releases", url("https://dl.bintray.com/content/sbt/sbt-plugin-releases"))(
  //  Resolver.ivyStylePatterns
  //)
)

addSbtPlugin("com.typesafe.sbt" % "sbt-ghpages" % "0.6.3")

addSbtPlugin("com.eed3si9n" % "sbt-unidoc" % "0.4.3")

addSbtPlugin("io.crashbox" % "sbt-gpg" % "0.2.1")

addSbtPlugin("org.xerial.sbt" % "sbt-sonatype" % "3.9.2")

// scoverage and coveralls deps are at old versions to avoid a bug in the current versions
// update these when this fix is released:  https://github.com/scoverage/sbt-coveralls/issues/73
//addSbtPlugin("org.scoverage" % "sbt-scoverage" % "1.0.4")

//addSbtPlugin("org.scoverage" % "sbt-coveralls" % "1.0.0")

//addSbtPlugin("org.scalastyle" %% "scalastyle-sbt-plugin" % "0.6.0")
