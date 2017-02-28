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

    <!-- Adding jar -->
    <dependencySets>
        <dependencySet>
            <useProjectArtifact>false</useProjectArtifact>
            <useTransitiveDependencies>false</useTransitiveDependencies>
            <unpack>false</unpack>
            <includes>
                <include>${project.groupId}:jom-sources:jar:*</include>
            </includes>
            <outputDirectory/>
            <outputFileNameMapping>jom-${project.version}.jar</outputFileNameMapping>
        </dependencySet>
    </dependencySets>

    <moduleSets>
        <!-- Javadoc -->
        <moduleSet>
            <useAllReactorProjects>true</useAllReactorProjects>
            <includes>
                <include>${project.groupId}:jom-javadoc</include>
            </includes>
            <binaries>
                <unpack>true</unpack>
                <attachmentClassifier>javadoc</attachmentClassifier>
                <outputDirectory>javadoc</outputDirectory>
                <includeDependencies>false</includeDependencies>
            </binaries>
        </moduleSet>
        <!-- Adding sources -->
        <moduleSet>
            <useAllReactorProjects>true</useAllReactorProjects>
            <includes>
                <include>${project.groupId}.jom-sources</include>
            </includes>
            <sources>
                <includeModuleDirectory>false</includeModuleDirectory>
                <fileSets>
                    <fileSet>
                        <outputDirectory>src</outputDirectory>
                    </fileSet>
                </fileSets>
            </sources>
        </moduleSet>
    </moduleSets>
</assembly>