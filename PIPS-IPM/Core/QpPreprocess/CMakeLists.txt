add_library(ooqpqpscaler STATIC QpScaler.C)
add_library(ooqpequistochscaler STATIC EquiStochScaler.C GeoStochScaler.C QpScaler.C StochScaler.C)
add_library(ooqpqppresolver STATIC QpPresolver.C)
add_library(ooqppqpostsolver STATIC QpPostsolver.C )
#add_library(   StochPresolverSingletonColumns.C )
add_library(ooqpqpstochpresolver STATIC StochRowStorage.C StochColumnStorage.C QpPostsolver.C StochPostsolver.C StochPresolver.C QpPresolver.C StochPresolverBase.C
 StochPresolverModelCleanup.C PresolveData.C StochPresolverBoundStrengthening.C StochPresolverColumnFixation.C StochPresolverParallelRows.C
 StochPresolverSingletonRows.C StochPresolverSingletonColumns.C )
