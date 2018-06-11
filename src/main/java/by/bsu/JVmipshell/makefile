INS=install
RM=rm -f
JHOME=/usr/lib/jvm/java
MIP_DIR=$(MIPCLDIR)
MIPCL_JAVA_DIR=$(MIP_DIR)/mipcl-java
INSDIR=$(MIPCL_JAVA_DIR)
TARGET=mipcl.jar

JSRC_PATH=$(MIPCL_JAVA_DIR)/sources/by/bsu/JVmipshell
CLASS_PATH=$(MIPCL_JAVA_DIR)/classes
vpath %.java $(JSRC_PATH)
JSRC=Var.java VarVector.java LinSum.java Ctr.java Function.java LPshell.java MIPshell.java
JCLS=$(JSRC:.java=.class)
%.class : %.java
	javac -cp $(CLASS_PATH) $< -d $(CLASS_PATH)
#
build: $(JCLS)
	jar cf $(INSDIR)/$(TARGET) -C $(CLASS_PATH) by
# java dependencies
Var.class: Var.java
VarVector.class: Var.java VarVector.java
LinSum.class: Var.java LinSum.java
Ctr.class: Var.java LinSum.java Ctr.java
Function.class: Var.java Function.java
LPshell.class: Var.java VarVector.java LinSum.java Ctr.java Function.java LPshell.java
MIPshell.class: Var.java VarVector.java LinSum.java Ctr.java Function.java MIPshell.java
.PHONY: clean
clean:
	$(RM) $(JCLS)	
test:

