"""Minimal rockmoc quickstart example."""

from rockmoc import DesignSpec, PerfectGas, design_ideal, design_mln, design_tic


def main() -> None:
    gas = PerfectGas(gamma=1.25, R=360.0)

    mln_spec = DesignSpec(
        gas=gas,
        throat_radius=0.05,  # m
        exit_mach=3.0,
        n_characteristics=40,
    )
    mln = design_mln(mln_spec)
    print("MLN:", mln.summary())

    ideal_spec = DesignSpec(
        gas=gas,
        throat_radius=0.05,
        exit_mach=3.0,
        n_characteristics=40,
        throat_rc_over_rt=1.5,
    )
    ideal = design_ideal(ideal_spec)
    print("Ideal:", ideal.summary())

    tic_spec = DesignSpec(
        gas=gas,
        throat_radius=0.05,
        exit_mach=3.0,
        n_characteristics=40,
        truncate_fraction=0.8,
    )
    tic = design_tic(tic_spec)
    print("TIC:", tic.summary())

    mln.export_csv("example_mln_wall.csv")
    print("Wrote example_mln_wall.csv")


if __name__ == "__main__":
    main()
