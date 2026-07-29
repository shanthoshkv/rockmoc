"""Student-facing helpers: gases, reports, explainers, tables, limits."""

from rockmoc.student.explainers import (
    EXPLAINERS,
    all_explainers_markdown,
    explainer,
    explainer_markdown,
)
from rockmoc.student.gases import (
    air,
    describe_example_gases,
    example_gas,
    example_gas_info,
    list_example_gases,
    lox_ch4,
    lox_lh2,
    lox_rp1,
)
from rockmoc.student.limits import (
    KNOWN_LIMITS,
    LIMITS_SHORT,
    limits_block,
    limits_bullet_list,
)
from rockmoc.student.report import (
    design_report_dict,
    format_design_report,
    generate_lab_report,
)
from rockmoc.student.tables import (
    comparison_to_csv,
    comparison_to_latex,
    comparison_to_markdown,
    comparison_to_string,
)

__all__ = [
    "EXPLAINERS",
    "KNOWN_LIMITS",
    "LIMITS_SHORT",
    "air",
    "all_explainers_markdown",
    "comparison_to_csv",
    "comparison_to_latex",
    "comparison_to_markdown",
    "comparison_to_string",
    "describe_example_gases",
    "design_report_dict",
    "example_gas",
    "example_gas_info",
    "explainer",
    "explainer_markdown",
    "format_design_report",
    "generate_lab_report",
    "limits_block",
    "limits_bullet_list",
    "list_example_gases",
    "lox_ch4",
    "lox_lh2",
    "lox_rp1",
]
