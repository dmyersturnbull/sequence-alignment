organization := "com.github.dmyersturnbull"

name := "sequence-alignment"

version := "0.1"

javacOptions ++= Seq("-source", "1.8", "-target", "1.8", "-Xlint:all")

resolvers ++= Seq(
	"Sonatype OSS Snapshots" at "https://oss.sonatype.org/content/repositories/snapshots"
)

libraryDependencies ++= Seq(
	"junit" % "junit" % "4.12" % "test",
	"org.biojava" % "biojava-alignment" % "4.1.0",
	"org.slf4j" % "slf4j-api" % "1.7.12",
	"org.assertj" % "assertj-core" % "3.1.0" % "test",
	"org.assertj" % "assertj-guava" % "1.3.1" % " test",
	"org.mockito" % "mockito-core" % "2.0.26-beta" % "test",
	"com.google.guava" % "guava" % "18.0",
	"org.slf4j" % "slf4j-api" % "1.7.12",
	"com.google.code.findbugs" % "jsr305" % "3.0.0",
	"javax.validation" % "validation-api" % "1.1.0.Final"
)
