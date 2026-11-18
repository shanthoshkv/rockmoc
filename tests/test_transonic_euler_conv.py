from rockmoc.moc.transonic_euler import cuffel_wall, solve_throat


def test_euler_throat_converges_to_steady_state():
    # regression: the old static-p inflow BC went period-2 at CFL 0.5, residual stuck ~3e-2
    r = solve_throat(cuffel_wall(0.625, r_inlet=2.0), -1.8, 1.0, nx=60, nr=16,
                     max_iter=5000, tol=1e-10)
    assert r.residual < 1e-10 and r.iterations < 5000
