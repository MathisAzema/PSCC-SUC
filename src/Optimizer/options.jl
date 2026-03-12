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
    compute_radius_DRO
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
    compute_radius_DRO
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
    compute_radius_DRO
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
    compute_radius_DRO
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
    oracle_RO_DCA_problem,
    DCAlgo_RO,
    add_cut_RO,
    _add_optimality_cuts_RO,
    get_variables_RO,
    get_value_variables_RO,
    compute_uncertainty_RO,
    results_second_stage_SP,
    compute_radius_RO
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
    compute_radius_RO
)

Grb_l2_budget = Options(
    "Grb_l2_budget",
    master_DRO_l2_budget_problem,
    oracle_DRO_l2_budget_problem,
    second_stage_grb_l2_budget,
    add_cut_DRO,
    _add_optimality_cuts_DRO,
    get_variables_DRO_l2_budget,
    get_value_variables_DRO_l2_budget,
    compute_uncertainty_l2_budget,
    results_second_stage_grb_l2_budget,
    compute_radius_DRO
)

DCA_l2_budget = Options(
    "DCA_l2_budget",
    master_DRO_l2_budget_problem,
    oracle_DRO_budget_DCA_problem,
    DCAlgoDRObudget,
    add_cut_DRO_budget,
    _add_optimality_cuts_DRO_budget,
    get_variables_DRO_l2_budget,
    get_value_variables_DRO_l2_budget,
    compute_uncertainty_l2_budget,
    results_second_stage_dca_l2_budget,
    compute_radius_DRO
)

Moment_MILP = Options(
    "Moment_MILP",
    master_moment_problem,
    oracle_moment_problem,
    second_stage_moment,
    add_cut_moment,
    _add_optimality_cuts_moment,
    get_variables_moment,
    get_value_variables_moment,
    compute_uncertainty_moment,
    results_second_stage_SP,
    compute_radius_RO
)

Moment_DCA = Options(
    "Moment_DCA",
    master_moment_problem,
    oracle_moment_problem,
    DCAlgo_moment,
    add_cut_moment,
    _add_optimality_cuts_moment,
    get_variables_moment,
    get_value_variables_moment,
    compute_uncertainty_moment,
    results_second_stage_SP,
    compute_radius_RO
)