package org.dde.icedeflection.deflection.parameter;

/**
 * Psi class
 *
 * @author  Nikolai Osipov <nao99.dev@gmail.com>
 * @version 1.0.0
 * @since   2021-04-04
 */
public class Psi {
    /**
     * Number of psi
     */
    private final int number;

    /**
     * Lower boundary
     */
    private final double lowerBoundary;

    /**
     * Upper boundary
     */
    private final double upperBoundary;

    /**
     * Steps number
     */
    private final int stepsNumber;

    /**
     * Step
     */
    private final double step;

    /**
     * Psi constructor
     *
     * @param number        a number of psi
     * @param lowerBoundary a lower boundary
     * @param upperBoundary an upper boundary
     * @param stepsNumber   a steps number
     */
    public Psi(int number, double lowerBoundary, double upperBoundary, int stepsNumber) {
        this.number = number;
        this.lowerBoundary = lowerBoundary;
        this.upperBoundary = upperBoundary;
        this.stepsNumber = stepsNumber;
        this.step = (upperBoundary - lowerBoundary) / stepsNumber;
    }

    public int getNumber() {
        return number;
    }

    public double getLowerBoundary() {
        return lowerBoundary;
    }

    public double getUpperBoundary() {
        return upperBoundary;
    }

    public int getStepsNumber() {
        return stepsNumber;
    }

    public double getStep() {
        return step;
    }
}
