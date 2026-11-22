"""Annular aerospike example: full Lee/Vernacchia spike, then a 25 % truncated plug."""

from rockmoc import design_spike, design_spike_truncated


def main() -> None:
    spike = design_spike(
        Pc=5e6, Tc=3200.0, gamma=1.25, molar_mass=22.0,
        Re=0.15, expansion_ratio=25.0, Pa=0.0, n_points=5000,
    )
    print("Full spike:", spike.summary())
    trunc = design_spike_truncated(spike, length_fraction=0.25)
    print("Truncated (25 %):", trunc.summary())


if __name__ == "__main__":
    main()
