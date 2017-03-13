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

    <!-- Adding project info -->
    <fileSets>
        <fileSet>
            <directory>${project.parent.basedir}</directory>
            <outputDirectory/>
            <includes>
                <include>README*</include>
                <include>LICENSE*</include>
                <include>NOTICE*</include>
            </includes>
        </fileSet>
    </fileSets>

    <dependencySets>
        <!-- Adding jar -->
        <dependencySet>
            <useProjectArtifact>false</useProjectArtifact>
            <useTransitiveDependencies>false</useTransitiveDependencies>
            <useStrictFiltering>true</useStrictFiltering>
            <unpack>false</unpack>
            <includes>
                <include>${project.groupId}:com.jom-sources:*</include>
            </includes>
            <outputDirectory/>
            <outputFileNameMapping>com.jom-${project.version}.jar</outputFileNameMapping>
        </dependencySet>
    </dependencySets>

    <moduleSets>
        <!-- Javadoc -->
        <moduleSet>
            <useAllReactorProjects>true</useAllReactorProjects>
            <includes>
                <include>${project.groupId}:com.jom-javadoc</include>
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
                <include>${project.groupId}:com.jom-sources</include>
            </includes>
            <sources>
                <includeModuleDirectory>false</includeModuleDirectory>
                <fileSets>
                    <fileSet>
                        <directory>src/main/java</directory>
                        <includes>
                            <include>**</include>
                        </includes>
                        <outputDirectory>src</outputDirectory>
                    </fileSet>
                </fileSets>
            </sources>
        </moduleSet>
    </moduleSets>
</assembly>