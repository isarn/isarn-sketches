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

// xsbt clean unidoc previewSite
// xsbt clean unidoc ghpagesPushSite
// xsbt -Dsbt.global.base=/home/eje/.sbt/sonatype +isarn_sketches/publish
// publish isarn-sketches-java for exactly one scala version:
// xsbt -Dsbt.global.base=/home/eje/.sbt/sonatype isarn_sketches_java/publish

scalaVersion := "2.11.12"

crossScalaVersions := Seq("2.11.12", "2.12.6")

// these do not "inherit" when defined at top level, so
// define them here for inclusion in each subproject.
// This also worked: 'xxx in ThisProject := yyy', but you have to do it
// for each setting below, so this seemed a bit cleaner
def publishSettings = Seq(
  version := "0.1.3-SNAPSHOT-pr12",
  //isSnapshot := true,
  //publishConfiguration := publishConfiguration.value.withOverwrite(true),
  //publishLocalConfiguration := publishLocalConfiguration.value.withOverwrite(true),
  organization := "org.isarnproject",
  pomIncludeRepository := { _ => false },
  publishMavenStyle := true,
  publishTo := {
    val nexus = "https://oss.sonatype.org/"
    if (isSnapshot.value)
      Some("snapshots" at nexus + "content/repositories/snapshots")
    else
      Some("releases"  at nexus + "service/local/staging/deploy/maven2")
  },
  licenses += ("Apache-2.0", url("http://opensource.org/licenses/Apache-2.0")),
  homepage := Some(url("https://github.com/isarn/isarn-sketches")),
  scmInfo := Some(
    ScmInfo(
      url("https://github.com/isarn/isarn-sketches"),
      "scm:git@github.com:isarn/isarn-sketches.git"
    )
  ),
  developers := List(
    Developer(
      id    = "erikerlandson",
      name  = "Erik Erlandson",
      email = "eje@redhat.com",
      url   = url("https://erikerlandson.github.io/")
    )
  )
)

compileOrder := CompileOrder.JavaThenScala

javacOptions ++= Seq()

scalacOptions ++= Seq("-unchecked", "-deprecation", "-feature")

scalacOptions in (Compile, doc) ++= Seq("-doc-root-content", baseDirectory.value+"/root-doc.txt")

enablePlugins(ScalaUnidocPlugin, JavaUnidocPlugin, GhpagesPlugin)

git.remoteRepo := "git@github.com:isarn/isarn-sketches.git"

siteSubdirName in ScalaUnidoc := "scala/api"

siteSubdirName in JavaUnidoc := "java/api"

addMappingsToSiteDir(mappings in (ScalaUnidoc, packageDoc), siteSubdirName in ScalaUnidoc)

addMappingsToSiteDir(mappings in (JavaUnidoc, packageDoc), siteSubdirName in JavaUnidoc)

// tell unidoc to not do scala-doc for the isarn-sketches-java (javadoc will still get created)
unidocProjectFilter in (ScalaUnidoc, unidoc) := inAnyProject -- inProjects(isarn_sketches_java)

// this target needs to execute only once, at the top level
// turn it off for any sub-projects
def siteSubProjectSettings = Seq(
  previewSite := {}
)

// browser insisted on caching some older generated site at the default (4000)
previewFixedPort := Some(4444)

lazy val isarn_sketches_java = (project in file("isarn-sketches-java"))
  .settings(name := "isarn-sketches-java")
  .enablePlugins(GenJavadocPlugin, PublishJavadocPlugin)
  .settings(siteSubProjectSettings :_*)
  .settings(
    crossPaths := false,                            // drop off Scala suffix from artifact names
    autoScalaLibrary := false                       // exclude scala-library from dependencies
    )
  .settings(publishSettings :_*)

lazy val isarn_sketches = (project in file("."))
  .aggregate(isarn_sketches_java)
  .dependsOn(isarn_sketches_java)
  .settings(name := "isarn-sketches")
  .settings(
    // isarn_sketches_java needs to be published separately to work with 'crossPaths := false'
    aggregate in publish := false,
    libraryDependencies ++= Seq(
      "org.isarnproject" %% "isarn-algebra-api" % "0.0.3",
      "org.isarnproject" %% "isarn-collections" % "0.0.4",
      "org.isarnproject" %% "isarn-scalatest" % "0.0.3" % Test,
      "org.scalatest" %% "scalatest" % "3.0.5" % Test,
      "org.apache.commons" % "commons-math3" % "3.6.1" % Test)
      )
  .settings(publishSettings :_*)
