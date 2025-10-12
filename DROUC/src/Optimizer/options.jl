risk_neutral = Options(
    "risk_neutral",
    master_SP_problem,
    oracle_SP_problem,
    second_stage_SP,
    add_cut_SP,
    _add_optimality_cuts_SP,
    get_variables_SP,
    get_value_variables_SP,
    compute_uncertainty_SP,
    results_second_stage_SP,
    compute_radius_SP
)

Grb_l2 = Options(
    "Grb_l2",
    master_DRO_l2_problem,
    oracle_DRO_l2_problem,
    second_stage_grb_l2,
    add_cut_DRO,
    _add_optimality_cuts_DRO,
    get_variables_DRO_l2,
    get_value_variables_DRO_l2,
    compute_uncertainty_SP,
    results_second_stage_grb_l2,
    compute_radius_l2
)

DCA_l2 = Options(
    "DCA_l2",
    master_DRO_l2_problem,
    oracle_DRO_l2_problem,
    DCAlgo,
    add_cut_DRO,
    _add_optimality_cuts_DRO,
    get_variables_DRO_l2,
    get_value_variables_DRO_l2,
    compute_uncertainty_l2,
    results_second_stage_dca_l2,
    compute_radius_l2
)

Grb_l1 = Options(
    "Grb_l1",
    master_DRO_l1_problem,
    oracle_DRO_l1_problem,
    second_stage_grb_l1,
    add_cut_DRO,
    _add_optimality_cuts_DRO,
    get_variables_DRO_l1,
    get_value_variables_DRO_l1,
    compute_uncertainty_SP,
    results_second_stage_grb_l1,
    compute_radius_l1
)

AVAR = Options(
    "AVAR",
    master_AVAR_problem,
    oracle_SP_problem,
    second_stage_SP,
    add_cut_AVAR,
    _add_optimality_cuts_AVAR,
    get_variables_AVAR,
    get_value_variables_AVAR,
    compute_uncertainty_SP,
    results_second_stage_SP,
    compute_radius_AVAR
)

DCA_l1 = Options(
    "DCA_l1",
    master_DRO_l1_problem,
    oracle_DRO_l2_problem,
    DCAlgo,
    add_cut_DRO,
    _add_optimality_cuts_DRO,
    get_variables_DRO_l1,
    get_value_variables_DRO_l1,
    compute_uncertainty_l1,
    results_second_stage_dca_l1,
    compute_radius_l1
)

KL = Options(
    "KL",
    master_KL_problem,
    oracle_SP_problem,
    second_stage_SP,
    add_cut_KL,
    _add_optimality_cuts_KL,
    get_variables_KL,
    get_value_variables_KL,
    compute_uncertainty_SP,
    results_second_stage_SP,
    compute_radius_KL
)

RO_DCA = Options(
    "RO_DCA",
    master_RO_problem_benders,
    oracle_RO_problem,
    DCAlgo_RO,
    add_cut_RO,
    _add_optimality_cuts_RO,
    get_variables_RO,
    get_value_variables_RO,
    compute_uncertainty_RO,
    results_second_stage_SP,
    compute_radius_KL
)

RO_MILP = Options(
    "RO_MILP",
    master_RO_problem_benders,
    oracle_RO_problem,
    second_stage_RO,
    add_cut_RO,
    _add_optimality_cuts_RO,
    get_variables_RO,
    get_value_variables_RO,
    compute_uncertainty_RO,
    results_second_stage_SP,
    compute_radius_KL
)