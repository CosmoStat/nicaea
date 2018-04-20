import os
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies

mb = module_builder.module_builder_t( [r"cosmo.hpp", "../cosmo/cosmo.h"]
                                      , gccxml_path=r"" 
                                      , include_paths=['../cosmo', '../coyote', '../halomodel', '../tools', '../cosmo']
                                      , define_symbols=['_NOFFTW_'] )
mb.decls().exclude()
mb.classes('cosmo').include()
mb.classes('cosmology').include()

mb.build_code_creator(module_name='pynicaea')

mb.write_module('./bindings.cpp')
