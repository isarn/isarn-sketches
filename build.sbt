// xsbt clean unidoc previewSite
// xsbt clean unidoc ghpagesPushSite
// xsbt -Dsbt.global.base=/home/eje/.sbt/sonatype +publish

//name := "isarn-sketches"

organization := "org.isarnproject"

version := "0.1.3-SNAPSHOT"

scalaVersion := "2.11.12"

crossScalaVersions := Seq("2.11.12", "2.12.6")

pomIncludeRepository := { _ => false }

publishMavenStyle := true

publishTo := {
  val nexus = "https://oss.sonatype.org/"
  if (isSnapshot.value)
    Some("snapshots" at nexus + "content/repositories/snapshots")
  else
    Some("releases"  at nexus + "service/local/staging/deploy/maven2")
}

licenses += ("Apache-2.0", url("http://opensource.org/licenses/Apache-2.0"))

homepage := Some(url("https://github.com/isarn/isarn-sketches"))

scmInfo := Some(
  ScmInfo(
    url("https://github.com/isarn/isarn-sketches"),
    "scm:git@github.com:isarn/isarn-sketches.git"
  )
)

developers := List(
  Developer(
    id    = "erikerlandson",
    name  = "Erik Erlandson",
    email = "eje@redhat.com",
    url   = url("https://erikerlandson.github.io/")
  )
)

compileOrder := CompileOrder.JavaThenScala

javacOptions ++= Seq()

scalacOptions ++= Seq("-unchecked", "-deprecation", "-feature")

scalacOptions in (Compile, doc) ++= Seq("-doc-root-content", baseDirectory.value+"/root-doc.txt")

enablePlugins(ScalaUnidocPlugin, JavaUnidocPlugin, GenJavadocPlugin, PublishJavadocPlugin, GhpagesPlugin)

siteSubdirName in ScalaUnidoc := "latest/api"

siteSubdirName in JavaUnidoc := "java/api"

addMappingsToSiteDir(mappings in (ScalaUnidoc, packageDoc), siteSubdirName in ScalaUnidoc)

addMappingsToSiteDir(mappings in (JavaUnidoc, packageDoc), siteSubdirName in JavaUnidoc)

git.remoteRepo := "git@github.com:isarn/isarn-sketches.git"

lazy val isarn_sketches_java = (project in file("isarn-sketches-java"))
  .settings(name := "isarn-sketches-java")
  .settings(crossPaths := false) // drop off Scala suffix from artifact names
  .settings(autoScalaLibrary := false) // exclude scala-library from dependencies

lazy val isarn_sketches = (project in file("."))
  .aggregate(isarn_sketches_java)
  .dependsOn(isarn_sketches_java)
  .settings(name := "isarn-sketches")
  .settings(
    libraryDependencies ++= Seq(
      "org.isarnproject" %% "isarn-algebra-api" % "0.0.3",
      "org.isarnproject" %% "isarn-collections" % "0.0.4",
      "org.isarnproject" %% "isarn-scalatest" % "0.0.3" % Test,
      "org.scalatest" %% "scalatest" % "3.0.5" % Test,
      "org.apache.commons" % "commons-math3" % "3.6.1" % Test))
