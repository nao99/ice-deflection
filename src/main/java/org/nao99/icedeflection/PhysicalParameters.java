package org.nao99.icedeflection;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

/**
 * PhysicalParameters class <br>
 *
 * This class contains all initial physical parameters
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-02-28
 */
public class PhysicalParameters {
    private static final double G = 9.8;

    /**
     * Ice place's thickness
     */
    private final double hi;

    /**
     * Half width of a canal (2L is a common width)
     */
    private final double L;

    /**
     * Canal's depth
     */
    private final double H;

    /**
     * Young's modulus
     */
    private final double E;

    /**
     * Ice's viscosity
     */
    private final double nu;

    /**
     * Delay time
     */
    private final double tau;

    /**
     * Ice plate's hardness
     */
    private final double D;

    /**
     * Canal's density
     */
    private final double pi;

    /**
     * Liquid's pressure to a lower surface of an ice plate
     */
    private final double pl;

    /**
     * Moving speed of an ice place
     */
    private final double U;

    /**
     * Weight of ice per unit area
     */
    private final double M;

    /**
     * Canal's depth without any dimension
     */
    private final double h;

    /**
     * Dimensionless parameter
     */
    private final double beta;

    /**
     * Dimensionless parameter
     */
    private final double epsilon;

    /**
     * Dimensionless parameter
     */
    private final double alpha;

    /**
     * Froude's number
     */
    private final double Fr;

    /**
     * PhysicalParameters constructor
     *
     * @param hi  an ice place's thickness
     * @param L   a half width of a canal
     * @param H   a canal's depth
     * @param E   a Young's modulus
     * @param nu  an ice's viscosity
     * @param tau a delay time
     * @param pi  a canal's density
     * @param pl  a liquid's pressure to a lower surface of an ice plate
     * @param U   a moving speed of an ice place
     */
    public PhysicalParameters(
        double hi,
        double L,
        double H,
        double E,
        double nu,
        double tau,
        double pi,
        double pl,
        double U
    ) {
        this.hi = hi;
        this.L = L;
        this.H = H;
        this.E = E * pow(10, 9);
        this.D = (this.E * pow(hi, 3)) / (12 * (1 - pow(nu, 2)));
        this.nu = nu;
        this.tau = tau;
        this.pi = pi;
        this.pl = pl;
        this.U = U;

        this.M = pi * hi;
        this.h = H / L;
        this.beta = D / (pl * G * pow(L, 4));
        this.epsilon = (tau * U) / L;
        this.alpha = M / (pl * L);
        this.Fr = U / sqrt(G * H);
    }

    public double getHi() {
        return hi;
    }

    public double getL() {
        return L;
    }

    public double getH() {
        return H;
    }

    public double getBeta() {
        return beta;
    }

    public double getEpsilon() {
        return epsilon;
    }

    public double getAlpha() {
        return alpha;
    }

    public double getFr() {
        return Fr;
    }

    public double getE() {
        return E;
    }

    public double getNu() {
        return nu;
    }

    public double getTau() {
        return tau;
    }

    public double getD() {
        return D;
    }

    public double getPi() {
        return pi;
    }

    public double getPl() {
        return pl;
    }

    public double getU() {
        return U;
    }

    public double getM() {
        return M;
    }

    public double getHDimensionless() {
        return h;
    }
}
