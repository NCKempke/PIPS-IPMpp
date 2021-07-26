### possibly add other readers - on the pyomo branch there is the pyomo dir here as well ! don't remove this file :P
set(includes_for_interface
        ../Core/Preprocessing
        ../Core/Interface
        ../Core/InteriorPointMethod
        ../Core/Utilities
        ../Core/Base
        ../Core/Options
        )

add_subdirectory(CallbackExample)
add_subdirectory(gams/gmspips)
