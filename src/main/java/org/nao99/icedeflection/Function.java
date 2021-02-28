package org.nao99.icedeflection;

/**
 * Function interface
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-22
 */
public interface Function {
    /**
     * Solves the equation
     *
     * @param x an x param
     * @return a solution
     */
    double solve(double x);
}
