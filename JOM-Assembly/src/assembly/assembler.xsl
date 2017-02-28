<assembly xmlns="http://maven.apache.org/ASSEMBLY/2.0.0"
          xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
          xsi:schemaLocation="http://maven.apache.org/ASSEMBLY/2.0.0 http://maven.apache.org/xsd/assembly-2.0.0.xsd">

    <id>distribution</id>

    <formats>
        <format>zip</format>
        <format>tar.gz</format>
        <format>tar.bz2</format>
    </formats>

    <!-- Building ZIP package -->

    <!-- Adding project info -->
    <fileSets>
        <fileSet>
            <directory>${project.basedir}</directory>
            <outputDirectory/>
            <includes>
                <include>README*</include>
                <include>LICENSE*</include>
                <include>NOTICE*</include>
            </includes>
        </fileSet>
    </fileSets>

    <!-- Adding sources -->
    <dependencySets>
        <dependencySet>
            <useProjectArtifact>false</useProjectArtifact>
            <useTransitiveDependencies>false</useTransitiveDependencies>
            <unpack>false</unpack>
            <includes>
                <include>${project.groupId}:jom-sources:jar:*</include>
            </includes>
            <outputDirectory/>
            <outputFileNameMapping>defaultNetworkDesign.jar</outputFileNameMapping>
        </dependencySet>
    </dependencySets>

    <!--<fileSet>-->
    <!--<directory>${project.build.directory}</directory>-->
    <!--<outputDirectory></outputDirectory>-->
    <!--<includes>-->
    <!--<include>*.jar</include>-->
    <!--</includes>-->
    <!--</fileSet>-->
    <!--<fileSet>-->
    <!--<directory>${project.basedir}/src/main/java</directory>-->
    <!--<outputDirectory>src</outputDirectory>-->
    <!--<includes>-->
    <!--<include>**</include>-->
    <!--</includes>-->
    <!--</fileSet>-->
    <!--<fileSet>-->
    <!--<directory>${project.build.directory}/site/apidocs</directory>-->
    <!--<outputDirectory>javadoc</outputDirectory>-->
    <!--<includes>-->
    <!--<include>**</include>-->
    <!--</includes>-->
    <!--</fileSet>-->
    <!--</fileSets>-->
</assembly>