  curvature_undulation_coupling:
    calculation: 
      curvature_undulation_coupling_dex:
        design: v3
      import_readymade_meso_v1_nanogel: {}
      import_readymade_meso_v1_membrane: {}
    collections: all
    specs: 
      # alternate names to match the calculation inputs above
      calcname: curvature_undulation_coupling_dex
      protein_abstractor_name: import_readymade_meso_v1_nanogel
      undulations_name: import_readymade_meso_v1_membrane
      # alternate loaders for interpreting the nanogel positions
      loader: 
        module: codes.curvature_coupling_loader_dex
        function: prepare_mesoscale_postdat
      # which plots to make
      routine:
        # - curvature_field_review
        - individual_reviews
