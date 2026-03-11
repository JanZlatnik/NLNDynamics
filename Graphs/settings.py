import numpy as np


MODELS = {

    # -------------------------------------------------------------------------------
    '2DModel2' : {
        # --- Potential settigns ---
        'POT_x_range': [0.0, 10],
        'POT_y_range': [-5.0, 2.0],
        'POT_nv_V0': 16,
        'POT_nv_Vd': 0,
        'POT_nv_VLCP': 6,
        'POT_n_ryd': 15,


        # --- DR Plot Settings ---
        'DR_2D_path': r"./DATA/drve_cs_model2.txt",
        'DR_Rydberg_E_path': r"/work/zlatnikj/DiscreteStateRuns/2D_M2ADS_FIT/DATA/RydbergStates/eigE_PHP.txt",
        'DR_Rydberg_R_path': r"/work/zlatnikj/DiscreteStateRuns/2D_M2ADS_FIT/DATA/RydbergStates/Rvalues.txt",
        'DR_x_range': [0.0, 1.5],
        'DR_y_range': [1e-4, 1e4],

        'DR_x_range_close': [[0.0, 0.3], [0.65, 0.95]],
        'DR_y_range_close': [[1e-2, 1e3], [1e-3, 1e2]],
        'DR_legend_loc': ['lower right', 'lower right'],
        'DR_legend_series_loc': ['lower left', 'lower left'],
        'DR_legend_frame': [True, True],

        'DR_v_plotted': [list(range(1, 101)), list(range(1, 101))],
        'DR_n_max': 200,  
        'DR_annotate_n': 10,       
        'DR_n_start': 1,          
        'DR_legend_prefix': r'$\nu=', 
        'DR_colormap': 'viridis',



        # --- VE Plot Settings ---
        'VE_v_plot': [0, 1, 2],
        'VE_x_range_close': {
            0: [[0.0, 1.5]],         
            1: [[0.25, 0.6]],         
            2: [[0.48, 0.8]]       
        },
        'VE_y_range_close_linear': {
            0: [[-10.0, 300.0]],
            1: [[-0.1, 5.0]],
            2: [[-0.5, 25.0]]
        },
        'VE_y_range_close_log': {
            0: [[1e-2, 1e4]],
            1: [[10**(-1.5), 10**(1.5)]],
            2: [[1e-3, 1e2]]
        },
        'VE_v_plotted': {
            0: [list(range(1, 101))],  
            1: [list(range(1, 101))],     
            2: [list(range(1, 101))]         
        },
        
        'VE_legend_loc_lin': {
            0: ['upper right'], 1: ['upper left'], 2: ['upper left']
        },
        'VE_legend_loc_log': {
            0: ['upper right'], 1: ['upper left'], 2: ['lower left']
        },
        
        'VE_legend_series_loc_lin': {
            0: ['upper left'], 1: ['lower left'], 2: ['lower left']
        },
        'VE_legend_series_loc_log': {
            0: ['lower left'], 1: ['lower left'], 2: ['lower left']
        },
        'VE_legend_frame': {
            0: [True], 1: [True], 2: [True]
        },
        
        'VE_n_max': 200,
        'VE_annotate_n': 10,
        'VE_n_start': 1,
        'VE_legend_prefix': r'$\nu=',
        'VE_colormap': 'viridis',


    },


     # -------------------------------------------------------------------------------
    '2DModel1' : {
        # --- Potential settigns ---
        'POT_x_range': [0.0, 10],
        'POT_y_range': [-5.0, 2.0],
        'POT_nv_V0': 16,
        'POT_nv_Vd': 0,
        'POT_nv_VLCP': 6,
        'POT_n_ryd': 15,


        # --- DR Plot Settings ---
        'DR_2D_path': r"./DATA/drve_cs_model1.txt",
        'DR_Rydberg_E_path': r"/work/zlatnikj/DiscreteStateRuns/2D_M1ADS_FIT/DATA/RydbergStates/eigE_PHP.txt",
        'DR_Rydberg_R_path': r"/work/zlatnikj/DiscreteStateRuns/2D_M1ADS_FIT/DATA/RydbergStates/Rvalues.txt",
        'DR_x_range': [0.0, 1.5],
        'DR_y_range': [1e-4, 1e4],

        'DR_x_range_close': [[0.0, 0.3], [0.65, 0.95]],
        'DR_y_range_close': [[1e-2, 1e3], [1e-3, 1e2]],
        'DR_legend_loc': ['lower right', 'lower right'],
        'DR_legend_series_loc': ['lower left', 'lower left'],
        'DR_legend_frame': [True, True],

        'DR_v_plotted': [list(range(1, 101)), list(range(1, 101))],
        'DR_n_max': 200,  
        'DR_annotate_n': 10,       
        'DR_n_start': 1,          
        'DR_legend_prefix': r'$\nu=', 


    },

  
}