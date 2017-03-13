<assembly xmlns="http://maven.apache.org/ASSEMBLY/2.0.0"
          xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
          xsi:schemaLocation="http://maven.apache.org/ASSEMBLY/2.0.0 http://maven.apache.org/xsd/assembly-2.0.0.xsd">

    <id>distribution</id>

    <formats>
        <format>zip</format>
        <format>tar.gz</format>
        <format>tar.bz2</format>
    </formats>

    <includeBaseDirectory>false</includeBaseDirectory>

    <!-- Building package -->

    <fileSets>
        <!-- Adding project info -->
        <fileSet>
            <directory>${project.parent.basedir}</directory>
            <outputDirectory/>
            <includes>
                <include>README*</include>
                <include>LICENSE*</include>
                <include>NOTICE*</include>
            </includes>
        </fileSet>
        <!-- Adding compiled file -->
        <fileSet>
            <directory>${project.build.directory}</directory>
            <outputDirectory/>
            <includes>
                <include>*.jar</include>
            </includes>
        </fileSet>
        <!-- Adding sources -->
        <fileSet>
            <directory>${project.basedir}/src/main/java</directory>
            <outputDirectory>src</outputDirectory>
            <includes>
                <include>**</include>
            </includes>
        </fileSet>
        <!-- Adding Javadoc -->
        <fileSet>
            <directory>${project.build.directory}/site/apidocs</directory>
            <outputDirectory>javadoc</outputDirectory>
            <includes>
                <include>**</include>
            </includes>
        </fileSet>
    </fileSets>
</assembly>