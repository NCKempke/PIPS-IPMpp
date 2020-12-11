
add_library(ooqpstoch STATIC 
   sFactory.C
   sFactoryAug.C
   #   sFactoryAugSchurLeaf.C
   #   sFactoryAugComm2SchurLeaf.C
   sData.C
   sLinsys.C
   sLinsysRoot.C
   sLinsysRootAug.C
   #   sLinsysRootComm2.C
   #   sLinsysRootAugComm2.C
   sLinsysLeaf.C
   sLinsysLeafSchurSlv.C
   sVars.C
   StochMonitor.C
   StochOptions.C
   sResiduals.C
   sTree.C
   sLinsysRootBordered.C
   sTreeCallbacks.C
   sLinsysLeafMumps.C
)
