"""Minimal rockmoc quickstart example (student-friendly)."""

from rockmoc import (
    DesignSpec,
    compare_nozzle_designs,
    design_conical,
    design_ideal,
    design_mln,
    design_tic,
    example_gas,
)


def main() -> None:
    # Homework-style gas + spec (Rt in mm)
    gas = example_gas("lox_rp1")
    spec = DesignSpec.student(Me=3.0, Rt_mm=50, n=40, gas_name="lox_rp1")

    mln = design_mln(spec)
    print(mln.report())
    print()

    ideal = design_ideal(
        DesignSpec(gas=gas, throat_radius=0.05, exit_mach=3.0, n_characteristics=40)
    )
    tic = design_tic(
        DesignSpec(
            gas=gas,
            throat_radius=0.05,
            exit_mach=3.0,
            n_characteristics=40,
            truncate_fraction=0.8,
        )
    )
    cone = design_conical(spec, half_angle_deg=15.0)
    print("Ideal L/Rt:", ideal.metrics.length / 0.05)
    print("TIC   L/Rt:", tic.metrics.length / 0.05)
    print("Cone  L/Rt:", cone.metrics.length / 0.05)

    cmp = compare_nozzle_designs(spec, include_conical=True)
    print()
    print(cmp.summary_table())

    mln.export_csv("example_mln_wall.csv")
    mln.lab_report("example_lab_report.md", student_name="Demo", course="quickstart")
    print("Wrote example_mln_wall.csv and example_lab_report.md")


if __name__ == "__main__":
    main()
