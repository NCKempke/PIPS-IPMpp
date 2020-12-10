add_library(ooqpequistochscaler STATIC EquiStochScaler.C GeoStochScaler.C StochScaler.C QpScaler.C )

add_library(ooqpqpstochpresolver STATIC StochRowStorage.C StochColumnStorage.C QpPostsolver.C StochPostsolver.C StochPresolver.C QpPresolver.C StochPresolverBase.C
 StochPresolverModelCleanup.C PresolveData.C StochPresolverBoundStrengthening.C StochPresolverColumnFixation.C StochPresolverParallelRows.C
StochPresolverSingletonRows.C StochPresolverSingletonColumns.C QpPresolver.C QpPostsolver.C)
